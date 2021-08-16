import os
from utils.transitive import TransitiveClosure

PDB_SEQDUMP_PATH = "data/SEQDUMP"
csv_path = "data/Ssym/Ssym_SCONES.csv"

RESSEQ_DUMP_PATH = "data/RESSEQ"

import pandas as pd

csv_data = pd.read_csv("data/Ssym/Ssym_SCONES.csv")
data = csv_data.T.to_dict()
data = { sample["id"] : sample for sample in data.values()}

transitive_closure = TransitiveClosure()
transitive_closure.load_dataset(data, RESSEQ_DUMP_PATH)
transitive_closure.transitive_closure()

def mergefn(sample1, sample2):
    sample3 = {
        "id" : sample1["id"] + '_' + sample2["id"],
        "pdb_id" : sample1["pdb_id"],
        "chain_id" : sample1["chain_id"],
        "uniprot_id" : sample1["uniprot_id"],
        "ref_residue" : sample1["ref_residue"],
        "seq_position" : sample1["seq_position"],
        "resseq_position" : sample1["resseq_position"],
        "mutated_residue" : sample2["mutated_residue"],
        "pH" : (sample1["pH"] + sample2["pH"])/2,
        "T" : (sample1["T"] + sample2["T"])/2,
        "ddG" : sample1["ddG"] + sample2["ddG"],
        "sequence" : sample1["sequence"]
    }

    from datasets.SsymDatasetGenerator import SsymDatasetGenerator
    assert(set(sample3.keys()) == set(SsymDatasetGenerator.SSYM_COLUMNS_LIST))
    return sample3

transitive_closure.export("data/Stransitive/Stransitive.csv", "data/Stransitive/Stransitive_idx.csv", mergefn)

import pandas as pd
csv_data = pd.read_csv("data/Stransitive/Stransitive.csv")
dataset = csv_data.T.to_dict()
dataset = { sample["id"] : sample for sample in dataset.values()}

from datasets.DatasetGenerator import DatasetGenerator
generator = DatasetGenerator()
generator.data = dataset
generator.validate_data(RESSEQ_DUMP_PATH)