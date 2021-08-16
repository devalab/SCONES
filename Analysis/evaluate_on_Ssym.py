import sys
import argparse
import csv

import numpy as np
import pandas as pd

def evaluate_on_Ssym(dataset: dict, symmetric_pairs: list, results: dict, fieldname, sign_flip):
    ddG_fwd_l, ddG_rev_l = [], []
    ddG_fwd_actual_l, ddG_rev_actual_l = [], []
    for A, B in symmetric_pairs:
        try:
            ddG_fwd = results[A][fieldname]
            ddG_rev = results[B][fieldname]

            ddG_fwd_l.append(ddG_fwd)
            ddG_fwd_actual_l.append(dataset[A]["ddG"])
            ddG_rev_l.append(ddG_rev)
            ddG_rev_actual_l.append(dataset[B]["ddG"])
        except Exception as ex:
            print("[FAILED] %s ~ %s" % (A, B))
            print("\t", type(ex).__name__, ex)

    ddG_fwd_l, ddG_fwd_actual_l = np.asarray(ddG_fwd_l), np.asarray(ddG_fwd_actual_l)
    ddG_rev_l, ddG_rev_actual_l = np.asarray(ddG_rev_l), np.asarray(ddG_rev_actual_l)

    if sign_flip:
        ddG_fwd_l = -ddG_fwd_l
        ddG_rev_l = -ddG_rev_l

    num_samples = ddG_rev_l.size
    Fc = np.corrcoef(ddG_fwd_l, ddG_fwd_actual_l)[0, 1]
    Rc = np.corrcoef(ddG_rev_l, ddG_rev_actual_l)[0, 1]
    Frmse = np.sqrt(np.mean((ddG_fwd_l - ddG_fwd_actual_l)**2))
    Rrmse = np.sqrt(np.mean((ddG_rev_l - ddG_rev_actual_l)**2))
    return Fc, Rc, Frmse, Rrmse, num_samples

def main():
    parser = argparse.ArgumentParser(description='Ssym Evaluator')
    parser.add_argument('--dataset_path', dest='dataset_path', type=str, help='path to Ssym dataset in CSV', required=True, default=None)
    parser.add_argument('--symmetric_pairs_list', dest='symmetric_pairs_list', type=str, help='path to symmetric pairs list in CSV', required=True, default=None)
    parser.add_argument('--results_path', dest='results_path', type=str, help='path to results CSV', required=True, default=None)
    parser.add_argument('--results_fieldname', dest='results_fieldname', type=str, help='results column that contains ddG values', default="ddG")
    parser.add_argument('--sign_flip', dest='sign_flip', type=bool, help='reverse ddG signs for evaluation', default=False)

    args = parser.parse_args()
    if args.dataset_path is None or args.symmetric_pairs_list is None or args.results_path is None:
        parser.print_help()
        sys.exit()

    dataset_path = args.dataset_path
    symmetric_pairs_list = args.symmetric_pairs_list
    results_path = args.results_path
    results_fieldname = args.results_fieldname
    sign_flip = args.sign_flip

    ssym = pd.read_csv(dataset_path)
    ssym = ssym.T.to_dict()
    ssym = { sample["id"] : sample for sample in ssym.values() }

    with open(symmetric_pairs_list, newline='') as f:
        reader = csv.reader(f)
        symmetric_pairs = list(reader)

    results = pd.read_csv(results_path)
    results = results.T.to_dict()
    results = { sample["id"] : sample for sample in results.values() }

    Fc, Rc, Frmse, Rrmse, N = evaluate_on_Ssym(ssym, symmetric_pairs, results, results_fieldname, sign_flip)
    print("Symmetric pairs considered: %d/%d" % (N, len(ssym) // 2))
    print("Forward correlation metric:", Fc)
    print("Reverse correlation metric:", Rc)
    print("Forward RMSE metric:", Frmse)
    print("Reverse RMSE metric:", Rrmse)


if __name__ == "__main__":
    main()
