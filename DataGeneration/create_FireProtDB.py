import numpy as np

def read_dataset(path: str) -> dict:
    from pandas import read_csv
    csv_data = read_csv(path)
    data = csv_data.T.to_dict()
    data = { sample["id"] : sample for sample in data.values() }
    return data

def merge_duplicates(source_path: str, target_path: str, ddG_stddev_cutoff: float=0.5, dTm_stddev_cutoff: float=1.0) -> None:
    print("MergeDuplicates(%s, %s)" % (source_path, target_path))

    def mergefn(group):
        ddG = [sample["ddG"] for sample in group if not np.isnan(sample["ddG"])]
        dTm = [sample["dTm"] for sample in group if not np.isnan(sample["dTm"])]
        pH = [sample["pH"] for sample in group if not np.isnan(sample["pH"])]

        sequence = group[0]["sequence"]
        seq_position = group[0]["seq_position"]
        mutated_residue = group[0]["mutated_residue"]
        ref_residue = group[0]["ref_residue"]
        pdb_id = group[0]["pdb_id"]
        chain = group[0]["chain_id"]

        for sample in group:
            assert(sequence == sample["sequence"])
            assert(seq_position == sample["seq_position"])
            assert(mutated_residue == sample["mutated_residue"])
            assert(ref_residue == sample["ref_residue"])
            assert(pdb_id == sample["pdb_id"])
            assert(chain == sample["chain_id"])

        if len(ddG) > 0 and np.std(ddG) > ddG_stddev_cutoff:
            ddG = []

        if len(dTm) > 0 and np.std(dTm) > dTm_stddev_cutoff:
            dTm = []

        if len(ddG) == 0 and len(dTm) == 0:
            return None

        merged_id = "_".join([sample["id"] for sample in group])
        merged_sample = group[0]
        merged_sample["id"] = merged_id
        merged_sample["pH"] = np.nan if len(pH) == 0 else np.mean(pH)
        merged_sample["ddG"] = np.nan if len(ddG) == 0 else np.mean(ddG)
        merged_sample["dTm"] = np.nan if len(dTm) == 0 else np.mean(dTm)
        return merged_sample

    data = read_dataset(source_path)

    from utils.merge_duplicates import MergeDuplicates
    merger = MergeDuplicates()
    merger.merge(data, RESSEQ_DUMP_PATH)
    merger.export(target_path, mergefn, columns=FireProtDBDatasetGenerator.FIREPROTDB_COLUMNS_LIST)

def symmetric_closure(source_path: str, target_path: str) -> None:
    print("SymmetricClosure(%s, %s)" % (source_path, target_path))

    def invertfn(sample):
        import copy
        fwd_sample = copy.deepcopy(sample)
        rev_sample = copy.deepcopy(sample)
        fwd_sample["id"] += "_F"
        rev_sample["id"] += "_R"

        ref_residue = sample["ref_residue"]
        mutated_residue = sample["mutated_residue"]
        seq_position = sample["seq_position"]

        rev_sequence = list(sample["sequence"])
        rev_sequence[seq_position] = mutated_residue
        rev_sequence = "".join(rev_sequence)

        rev_sample["sequence"] = rev_sequence
        rev_sample["ref_residue"] = mutated_residue
        rev_sample["mutated_residue"] = ref_residue
        rev_sample["ddG"] = -sample["ddG"]
        rev_sample["dTm"] = -sample["dTm"]
        return fwd_sample, rev_sample

    data = read_dataset(source_path)

    from utils.symmetric import SymmetricClosure
    symclosure = SymmetricClosure()
    symclosure.augment(data, invertfn)
    symclosure.export(target_path, columns=FireProtDBDatasetGenerator.FIREPROTDB_COLUMNS_LIST)

def transitive_closure(source_path: str, target_path_csv: str, target_path_idx: str) -> None:
    print("TransitiveClosure(%s, %s, %s)" % (source_path, target_path_csv, target_path_idx))

    def transitive_mergefn(sample1, sample2):
        import copy
        sample3 = copy.deepcopy(sample1)
        sample3["id"] = sample1["id"] + "_+_" + sample2["id"]
        sample3["mutated_residue"] = sample2["mutated_residue"]
        sample3["pH"] = (sample1["pH"] + sample2["pH"])/2
        sample3["ddG"] = sample1["ddG"] + sample2["ddG"]
        return sample3

    data = read_dataset(source_path)

    from utils.transitive import TransitiveClosure
    transclosure = TransitiveClosure()
    transclosure.load_dataset(data, RESSEQ_DUMP_PATH)
    transclosure.transitive_closure()
    transclosure.export(target_path_csv, target_path_idx, transitive_mergefn, columns=FireProtDBDatasetGenerator.FIREPROTDB_COLUMNS_LIST)

FireProtDB_CSV_PATH = "data/FireProtDB/fireprotdb_results.csv"
PDB_DUMP_PATH = "data/PDB"
LEGACY_PDB_DUMP_PATH = "data/LegacyPDB"
UNIPROT_DUMP_PATH = "data/UniProt"
RESSEQ_DUMP_PATH = "data/RESSEQ"

from datasets.FireProtDBDatasetGenerator import FireProtDBDatasetGenerator
generator = FireProtDBDatasetGenerator()
generator.read_FireProtDB(FireProtDB_CSV_PATH)
generator.download_pdb_files(PDB_DUMP_PATH)
generator.download_pdb_files(LEGACY_PDB_DUMP_PATH, file_format="pdb")
generator.download_uniprot(UNIPROT_DUMP_PATH)
generator.create_resseq(PDB_DUMP_PATH, RESSEQ_DUMP_PATH)
generator.calculate_resseq_from_seq(RESSEQ_DUMP_PATH)
generator.validate_data(RESSEQ_DUMP_PATH)
generator.export("data/FireProtDB/fireprotdb_cleaned.csv")

# merge_duplicates("data/FireProtDB/fireprotdb_cleaned.csv", "data/FireProtDB/fireprotdb_cleaned_unique.csv")
# symmetric_closure("data/FireProtDB/fireprotdb_cleaned_unique.csv", "data/FireProtDB/fireprotdb_cleaned_unique_symmetric.csv")

# transitive_closure("data/FireProtDB/fireprotdb_cleaned_unique.csv", "data/FireProtDB/fireprotdb_cleaned_unique_transitive.csv", "data/FireProtDB/fireprotdb_cleaned_unique_transitive_idx.csv")
# #transitive_closure("data/FireProtDB/fireprotdb_cleaned_symmetric_unique.csv", "data/FireProtDB/fireprotdb_cleaned_symmetric_unique_transitive.csv", "data/FireProtDB/fireprotdb_cleaned_symmetric_unique_transitive_idx.csv")
# transitive_closure("data/FireProtDB/fireprotdb_cleaned_unique_symmetric.csv", "data/FireProtDB/fireprotdb_cleaned_unique_symmetric_transitive.csv", "data/FireProtDB/fireprotdb_cleaned_unique_symmetric_transitive_idx.csv")
