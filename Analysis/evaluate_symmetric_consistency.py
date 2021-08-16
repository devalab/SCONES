import sys
import argparse
import csv

import numpy as np
import pandas as pd

def evaluate_symmetric_consistency(symmetric_pairs: list, results: dict, fieldname):
    ddG_XY, ddG_YX = [], []
    for A, B in symmetric_pairs:
        try:
            sampleA, sampleB = results[A], results[B]
            ddG_A, ddG_B = sampleA[fieldname], sampleB[fieldname]
            ddG_XY.append(ddG_A)
            ddG_YX.append(ddG_B)
        except Exception as ex:
            print("[FAILED] %s ~ %s" % (A, B))
            print("\t", type(ex).__name__, ex)

    ddG_XY = np.asarray(ddG_XY)
    ddG_YX = np.asarray(ddG_YX)

    num_samples = ddG_XY.size
    Rs = np.corrcoef(ddG_XY, ddG_YX)[0][1]
    ds = np.mean(ddG_XY + ddG_YX)
    Ds = np.sqrt(np.mean((ddG_XY + ddG_YX)**2))
    return Rs, ds, Ds, num_samples

def main():
    parser = argparse.ArgumentParser(description='Symmetric Consistency Evaluator')
    parser.add_argument('--symmetric_pairs_list', dest='symmetric_pairs_list', type=str, help='path to symmetric pairs list in CSV', required=True, default=None)
    parser.add_argument('--results_path', dest='results_path', type=str, help='path to results CSV', required=True, default=None)
    parser.add_argument('--results_fieldname', dest='results_fieldname', type=str, help='results column that contains ddG values', default="ddG")

    args = parser.parse_args()
    if args.symmetric_pairs_list is None or args.results_path is None:
        parser.print_help()
        sys.exit()

    symmetric_pairs_list = args.symmetric_pairs_list
    results_path = args.results_path
    results_fieldname = args.results_fieldname

    with open(symmetric_pairs_list, newline='') as f:
        reader = csv.reader(f)
        symmetric_pairs = list(reader)

    results = pd.read_csv(results_path)
    results = results.T.to_dict()
    results = { sample["id"] : sample for sample in results.values() }

    Rs, ds, Ds, N = evaluate_symmetric_consistency(symmetric_pairs, results, results_fieldname)
    print("Symmetric pairs considered: %d/%d" % (N, len(symmetric_pairs)))
    print("Rs correlation metric:", Rs)
    print("ds bias metric:", ds)
    print("Ds norm metric:", Ds)

if __name__ == "__main__":
    main()
