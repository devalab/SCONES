import os
import sys
import argparse
import csv

import numpy as np
import pandas as pd
import pymol
import seaborn as sns
import matplotlib.pyplot as plt

def analyze_structural_changes(dataset: dict, symmetric_pairs: list, pdb_fieldname: str, pdb_dump_path: str):
    ca_rmsd = []
    for A, B in symmetric_pairs:
        try:
            A_path = os.path.join(pdb_dump_path, dataset[A][pdb_fieldname] + ".cif")
            B_path = os.path.join(pdb_dump_path, dataset[B][pdb_fieldname] + ".cif")

            pymol.cmd.reinitialize()
            pymol.cmd.load(A_path, 'apo')
            pymol.cmd.load(B_path, 'holo')

            rmsd = pymol.cmd.align('holo' + '////CA', 'apo' +  '////CA', cycles=0, transform=0, object='aln')[3]
            ca_rmsd.append(rmsd)
        except Exception as ex:
            print("[FAILED] %s ~ %s" % (A, B))
            print("\t", type(ex).__name__, ex)

    ca_rmsd = np.asarray(ca_rmsd)

    mean = np.mean(ca_rmsd)
    median = np.median(ca_rmsd)
    print("mean:", mean)
    print("median:", median)
    print("more than one:", np.sum(ca_rmsd > 1))
    print("more than two:", np.sum(ca_rmsd > 2))
    print("more than five:", np.sum(ca_rmsd > 5))

    print(ca_rmsd)
    sns.histplot(ca_rmsd, bins=200, binrange=(0, np.max(ca_rmsd)))
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Analyzes structural changes between reference and mutant structure')
    parser.add_argument('--dataset_path', dest='dataset_path', type=str, help='path to the dataset with structures', required=True, default=None)
    parser.add_argument('--symmetric_pairs_list', dest='symmetric_pairs_list', type=str, help='path to symmetric pairs list in CSV', required=True, default=None)
    parser.add_argument('--pdb_fieldname', dest='pdb_fieldname', type=str, help='column that contains PDB ids', default="pdb_id")
    parser.add_argument('--pdb_dump_path', dest='pdb_dump_path', type=str, help='path to directory containing pdb ids', required=True)

    args = parser.parse_args()
    if args.dataset_path is None or args.symmetric_pairs_list is None:
        parser.print_help()
        sys.exit()

    dataset_path = args.dataset_path
    symmetric_pairs_list = args.symmetric_pairs_list
    pdb_fieldname = args.pdb_fieldname
    pdb_dump_path = args.pdb_dump_path

    dataset = pd.read_csv(dataset_path)
    dataset = dataset.T.to_dict()
    dataset = { sample["id"] : sample for sample in dataset.values() }

    with open(symmetric_pairs_list, newline='') as f:
        reader = csv.reader(f)
        symmetric_pairs = list(reader)

    analyze_structural_changes(dataset, symmetric_pairs, pdb_fieldname, pdb_dump_path)

if __name__ == "__main__":
    main()
