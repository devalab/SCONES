import sys
import argparse

import numpy as np
import pandas as pd

def evaluate_predictions(dataset: dict, results: dict, fieldname, sign_flip, min_magnitude_cutoff=None, max_magintude_cutoff=None):
    ddG_predicted, ddG_target = [], []
    for key, sample in dataset.items():
        if min_magnitude_cutoff and np.abs(sample["ddG"]) < min_magnitude_cutoff:
            continue
        if max_magintude_cutoff and np.abs(sample["ddG"]) > max_magintude_cutoff:
            continue
        try:
            ddG = results[key][fieldname]
            if sign_flip:
                ddG *= -1
            ddG_predicted.append(ddG)
            ddG_target.append(sample["ddG"])
        except Exception as ex:
            print("[FAILED] %s" % key)
            print("\t", type(ex).__name__, ex)

    ddG_predicted = np.asarray(ddG_predicted)
    ddG_target = np.asarray(ddG_target)

    num_samples = ddG_predicted.size
    pcc = np.corrcoef(ddG_predicted, ddG_target)[0, 1]
    mae = np.mean(np.abs(ddG_predicted - ddG_target))
    rmse = np.sqrt(np.mean((ddG_predicted - ddG_target)**2))
    return pcc, mae, rmse, num_samples

def main():
    parser = argparse.ArgumentParser(description='Basic Performance Evaluation')
    parser.add_argument('--dataset_path', dest='dataset_path', type=str, help='path to dataset in CSV', required=True, default=None)
    parser.add_argument('--results_path', dest='results_path', type=str, help='path to results CSV', required=True, default=None)
    parser.add_argument('--results_fieldname', dest='results_fieldname', type=str, help='results column that contains ddG values', default="ddG")
    parser.add_argument('--sign_flip', dest='sign_flip', type=bool, help='reverse ddG signs for evaluation', default=False)
    parser.add_argument('--min_magnitude_cutoff', dest='min_magnitude_cutoff', type=float, help='minimum experimental ddG magnitude for inclusion', default=None)
    parser.add_argument('--max_magnitude_cutoff', dest='max_magnitude_cutoff', type=float, help='maximum experimental ddG magnitude for inclusion', default=None)

    args = parser.parse_args()
    if args.dataset_path is None or args.results_path is None:
        parser.print_help()
        sys.exit()

    dataset_path = args.dataset_path
    results_path = args.results_path
    results_fieldname = args.results_fieldname
    sign_flip = args.sign_flip
    min_magnitude_cutoff = args.min_magnitude_cutoff
    max_magnitude_cutoff = args.max_magnitude_cutoff

    dataset = pd.read_csv(dataset_path, dtype = { "id" : str })
    dataset = dataset.T.to_dict()
    dataset = { sample["id"] : sample for sample in dataset.values() }

    results = pd.read_csv(results_path, dtype = { "id" : str })
    results = results.T.to_dict()
    results = { sample["id"] : sample for sample in results.values() }

    pcc, mae,  rmse, N = evaluate_predictions(dataset, results, results_fieldname, sign_flip, min_magnitude_cutoff, max_magnitude_cutoff)

    print("Overall Performance:")
    print("\tNumber of samples considered: %d/%d" % (N, len(dataset)))
    print("\tPCC:", pcc)
    print("\tMAE:", mae)
    print("\tRMSE:", rmse)
    print()

if __name__ == "__main__":
    main()
