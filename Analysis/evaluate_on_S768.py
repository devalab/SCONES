import sys
import argparse

import numpy as np
import pandas as pd

def evaluate_on_S768_categorywise(dataset: dict, results: dict, fieldname, sign_flip):
    def get_category(sample):
        return [cat for cat in sample["categories"].split('\'') if cat not in (' ', '')]

    category_names = set()
    for sample in dataset.values():
        category_names.update(get_category(sample))

    categorywise_data = {}
    for category in category_names:
        categorywise_data[category] = {
            "target" : [],
            "prediction" : []
        }

    for key, sample in dataset.items():
        try:
            ddG_predicted = results[key][fieldname]
            if sign_flip:
                ddG_predicted *= -1
            ddG_target = sample["ddG"]
            cats = get_category(sample)
            for cat in cats:
                categorywise_data[cat]["target"].append(ddG_target)
                categorywise_data[cat]["prediction"].append(ddG_predicted)
        except Exception as ex:
            print("[FAILED] %s" % key)
            print("\t", type(ex).__name__, ex)

    df = pd.DataFrame(0, index=category_names, columns=["num_samples", "pcc", "mae", "rmse"])
    df = df.astype({ "num_samples" : int, "pcc" : float, "mae" : float, "rmse" : float})
    for category in category_names:
        target_values = np.asarray(categorywise_data[category]["target"])
        predicted_values = np.asarray(categorywise_data[category]["prediction"])
        pcc = np.corrcoef(target_values, predicted_values)[0][1]
        mae = np.mean(np.abs(target_values - predicted_values))
        rmse = np.sqrt(np.mean((target_values - predicted_values)**2))
        df.at[category, "num_samples"] = len(target_values)
        df.at[category, "pcc"] = pcc
        df.at[category, "mae"] = mae
        df.at[category, "rmse"] = rmse

    df.sort_index(inplace=True)
    return df

def evaluate_on_S768(dataset: dict, results: dict, fieldname, sign_flip):
    ddG_predicted, ddG_target = [], []
    for key, sample in dataset.items():
        try:
            ddG = results[key][fieldname]
            if sign_flip:
                ddG *= -1
            ddG_predicted.append(ddG)
            ddG_target.append(sample["ddG"])
        except Exception as ex:
            print("[FAILED] %s" % key)
            print("\t", type(ex).__name__, ex)

    if len(ddG_predicted) == 0:
        raise Exception("No results were read. Please provide the correct fieldname for the column that contains the results.")

    ddG_predicted = np.asarray(ddG_predicted)
    ddG_target = np.asarray(ddG_target)

    num_samples = ddG_predicted.size
    pcc = np.corrcoef(ddG_predicted, ddG_target)[0, 1]
    rmse = np.sqrt(np.mean((ddG_predicted - ddG_target)**2))
    return pcc, rmse, num_samples

def main():
    parser = argparse.ArgumentParser(description='S768 Evaluator')
    parser.add_argument('--dataset_path', dest='dataset_path', type=str, help='path to S768 dataset in CSV', required=True, default=None)
    parser.add_argument('--results_path', dest='results_path', type=str, help='path to results CSV', required=True, default=None)
    parser.add_argument('--results_fieldname', dest='results_fieldname', type=str, help='results column that contains ddG values', default="ddG")
    parser.add_argument('--sign_flip', dest='sign_flip', type=bool, help='reverse ddG signs for evaluation', default=False)

    args = parser.parse_args()
    if args.dataset_path is None or args.results_path is None:
        parser.print_help()
        sys.exit()

    dataset_path = args.dataset_path
    results_path = args.results_path
    results_fieldname = args.results_fieldname
    sign_flip = args.sign_flip

    s768 = pd.read_csv(dataset_path)
    s768 = s768.T.to_dict()
    s768 = { int(sample["id"]) : sample for sample in s768.values() }

    results = pd.read_csv(results_path)
    results = results.T.to_dict()
    results = { int(sample["id"]) : sample for sample in results.values() }

    pcc, rmse, N = evaluate_on_S768(s768, results, results_fieldname, sign_flip)
    df = evaluate_on_S768_categorywise(s768, results, results_fieldname, sign_flip)

    print("Overall Performance:")
    print("\tNumber of samples considered: %d/%d" % (N, len(s768)))
    print("\tPCC:", pcc)
    print("\tRMSE:", rmse)
    print()

    print("Category-wise performance:")
    print(df)

if __name__ == "__main__":
    main()
