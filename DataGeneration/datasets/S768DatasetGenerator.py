from .DatasetGenerator import DatasetGenerator

class S768DatasetGenerator(DatasetGenerator):
    S768_COLUMNS_LIST=["id", "pdb_id", "chain_id", "ref_residue", "seq_position", "resseq_position", "mutated_residue", "ddG", "categories", "sequence"]
    
    def __init__(self):
        super().__init__()

    def read_S768(self, csv_path: str) -> None:
        from pandas import read_csv
        csv_data = read_csv(csv_path)
        for idx, sample in csv_data.iterrows():
            pdb_id, chain_id = sample["PDB ID"].lower(), sample["Chain"].upper()

            ref_residue = sample["Wild Type"]
            resseq_position = sample["Residue Number"]
            if resseq_position == "27D": # correction
                resseq_position = "27"
            resseq_position = int(resseq_position)
            mutated_residue = sample["Mutation"]

            self.add(idx,
                pdb_id=pdb_id,
                chain_id=chain_id,
                ref_residue=ref_residue,
                resseq_position=resseq_position,
                mutated_residue=mutated_residue,
                ddG=-float(sample["Experimental DDG"]),
                extra={ "categories" : sample["Classifiers"] }
            )

    def export(self, csv_path: str) -> None:
        super().export(csv_path, columns=self.S768_COLUMNS_LIST)
