import sys
import argparse
import csv

import numpy as np
import pandas as pd

def evaluate_transitive_consistency(transitive_pairs: list, results: dict, fieldname):
    ddG_XYYZ, ddG_XZ = [], []
    for A, B, C in transitive_pairs:
        try:
            sample_A, sample_B, sampleC = results[A], results[B], results[C]
            ddG_A, ddG_B, ddG_C = sample_A[fieldname], sample_B[fieldname], sampleC[fieldname]
            ddG_XYYZ.append(ddG_A + ddG_B)
            ddG_XZ.append(ddG_C)
        except Exception as ex:
            print("[FAILED] %s + %s -> %s" % (A, B, C))
            print("\t", type(ex).__name__, ex)

    ddG_XYYZ = np.asarray(ddG_XYYZ)
    ddG_XZ = np.asarray(ddG_XZ)

    num_samples = ddG_XZ.size
    Rt = np.corrcoef(ddG_XYYZ, ddG_XZ)[0, 1]
    Dt = np.sqrt(np.mean((ddG_XYYZ - ddG_XZ)**2))
    return Rt, Dt, num_samples

def main():
    parser = argparse.ArgumentParser(description='Transitive Consistency Evaluator')
    parser.add_argument('--transitive_tuples_list', dest='transitive_tuples_path', type=str, help='path to transitive tuples list in CSV', required=True, default=None)
    parser.add_argument('--results_path', dest='results_path', type=str, help='path to results CSV', required=True, default=None)
    parser.add_argument('--results_fieldname', dest='results_fieldname', type=str, help='results column that contains ddG values', default="ddG")

    args = parser.parse_args()
    if args.transitive_tuples_path is None or args.results_path is None:
        parser.print_help()
        sys.exit()

    transitive_tuples_path = args.transitive_tuples_path
    results_path = args.results_path
    results_fieldname = args.results_fieldname

    with open(transitive_tuples_path, newline='') as f:
        reader = csv.reader(f)
        transitive_pairs = list(reader)

    results = pd.read_csv(results_path)
    results = results.T.to_dict()
    results = { sample["id"] : sample for sample in results.values() }

    Rt, Dt, N = evaluate_transitive_consistency(transitive_pairs, results, results_fieldname)
    print("Transitive pairs considered: %d/%d" % (N, len(transitive_pairs)))
    print("Rt correlation metric:", Rt)
    print("Dt norm metric:", Dt)

if __name__ == "__main__":
    main()
