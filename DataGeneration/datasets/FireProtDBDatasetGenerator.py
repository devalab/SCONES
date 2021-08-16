from .DatasetGenerator import DatasetGenerator

class FireProtDBDatasetGenerator(DatasetGenerator):
    FIREPROTDB_COLUMNS_LIST=["id", "pdb_id", "chain_id", "ref_residue", "seq_position", "resseq_position", "mutated_residue", "pH", "ddG", "dTm", "sequence"]

    def __init__(self):
        super().__init__()

    def read_FireProtDB(self, csv_path: str, require_ddG: bool=False, require_dTm: bool=False) -> None:
        counts = {
            "residue_mismatch" : 0,
            "no_mutation" : 0,
            "ddG_not_available" : 0,
            "dTm_not_available" : 0,
            "no_ddG_or_dTm" : 0,
            "no_structure" : 0,
            "duplicate_sample" : 0
        }

        import numpy as np
        from pandas import read_csv
        csv_data = read_csv(csv_path)
        for _, sample in csv_data.iterrows():
            pdb_id = sample["pdb_id"]
            if not isinstance(pdb_id, str):
                counts["no_structure"] += 1
                continue

            sequence = sample["sequence"]
            seq_position = int(sample["position"]) - 1
            mutated_residue = sample["mutation"]
            ref_residue = sample["wild_type"]
            
            if sequence[seq_position] != ref_residue:
                counts["residue_mismatch"] += 1
                continue

            if ref_residue == mutated_residue:
                counts["no_mutation"] += 1
                continue

            if require_ddG and np.isnan(sample["ddG"]):
                counts["ddG_not_available"] += 1
                continue

            if require_dTm and np.isnan(sample["dTm"]):
                counts["dTm_not_available"] += 1
                continue

            if np.isnan(sample["ddG"]) and np.isnan(sample["dTm"]):
                counts["no_ddG_or_dTm"] += 1 
                continue

            try:
                self.add(sample["experiment_id"],
                    pdb_id=pdb_id,
                    chain_id=sample["chain"],
                    ref_residue=ref_residue,
                    sequence=sequence,
                    seq_position=seq_position,
                    mutated_residue=mutated_residue,
                    pH=sample["pH"],
                    ddG=-float(sample["ddG"]),
                    dTm=-float(sample["dTm"])
                )
            except ValueError as e:
                counts["duplicate_sample"] += 1 
                print(e)

        print(counts)

    def export(self, csv_path: str) -> None:
        super().export(csv_path, columns=self.FIREPROTDB_COLUMNS_LIST)
