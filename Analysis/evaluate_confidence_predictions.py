import sys
import argparse

import numpy as np
import pandas as pd

def evaluate_confidence_predictions(dataset: dict, results: dict, ddG_fieldname, sign_flip, conf_fieldname):
    ddG_predicted, ddG_predicted_conf, ddG_target = [], [], []
    for key, sample in dataset.items():
        try:
            ddG = results[key][ddG_fieldname]
            conf = results[key][conf_fieldname]
            if sign_flip:
                ddG *= -1
            ddG_predicted.append(ddG)
            ddG_predicted_conf.append(conf)
            ddG_target.append(sample["ddG"])
        except Exception as ex:
            print("[FAILED] %s" % key)
            print("\t", type(ex).__name__, ex)

    ddG_predicted = np.asarray(ddG_predicted)
    ddG_predicted_conf = np.asarray(ddG_predicted_conf)
    ddG_target = np.asarray(ddG_target)

    num_samples = ddG_predicted.size
    abs_error = np.abs(ddG_predicted - ddG_target)/(1e-5 + np.abs(ddG_target))
    pcc = np.corrcoef(abs_error, ddG_predicted_conf)[0, 1]
    return pcc, num_samples

def main():
    parser = argparse.ArgumentParser(description='Confidence Evaluation')
    parser.add_argument('--dataset_path', dest='dataset_path', type=str, help='path to dataset in CSV', required=True, default=None)
    parser.add_argument('--results_path', dest='results_path', type=str, help='path to results CSV', required=True, default=None)
    parser.add_argument('--ddG_fieldname', dest='ddG_fieldname', type=str, help='results column that contains ddG values', default="ddG")
    parser.add_argument('--conf_fieldname', dest='conf_fieldname', type=str, help='results column that contains ddG confidence values', default="conf")
    parser.add_argument('--sign_flip', dest='sign_flip', type=bool, help='reverse ddG signs for evaluation', default=False)

    args = parser.parse_args()
    if args.dataset_path is None or args.results_path is None:
        parser.print_help()
        sys.exit()

    dataset_path = args.dataset_path
    results_path = args.results_path
    ddG_fieldname = args.ddG_fieldname
    conf_fieldname = args.conf_fieldname
    sign_flip = args.sign_flip

    dataset = pd.read_csv(dataset_path, dtype = { "id" : str })
    dataset = dataset.T.to_dict()
    dataset = { sample["id"] : sample for sample in dataset.values() }

    results = pd.read_csv(results_path, dtype = { "id" : str })
    results = results.T.to_dict()
    results = { sample["id"] : sample for sample in results.values() }

    pcc, N = evaluate_confidence_predictions(dataset, results, ddG_fieldname, sign_flip, conf_fieldname)

    print("Overall Performance:")
    print("\tNumber of samples considered: %d/%d" % (N, len(dataset)))
    print("\tPCC:", pcc)
    print()

if __name__ == "__main__":
    main()
