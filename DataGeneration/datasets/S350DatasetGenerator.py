from .DatasetGenerator import DatasetGenerator

class S350DatasetGenerator(DatasetGenerator):
    S350_COLUMNS_LIST=["id", "pdb_id", "chain_id", "uniprot_id", "is_in_309", "ref_residue", "seq_position", "resseq_position", "mutated_residue", "pH", "T", "ddG", "sequence"]
    
    def __init__(self):
        super().__init__()

    def read_S350(self, csv_path: str) -> None:
        from pandas import read_csv
        csv_data = read_csv(csv_path)
        for idx, sample in csv_data.iterrows():
            pdb_id, chain_id = sample["PDB_CHAIN"][:4].lower(), sample["PDB_CHAIN"][4].upper()
            if pdb_id == "1rtp": # apply correction
                chain_id = '1'
            ref_residue = sample["WILD_RES"]
            resseq_position = int(sample["POSITION"])
            mutated_residue = sample["MUTANT_RES"]

            self.add(idx,
                pdb_id=pdb_id,
                chain_id=chain_id,
                ref_residue=ref_residue,
                resseq_position=resseq_position,
                mutated_residue=mutated_residue,
                pH=float(sample["PH"]),
                T=float(sample["TEMPERATURE"]),
                ddG=float(sample["EXP_DDG"])
            )

    def read_S350_PremPS(self, csv_path: str) -> None:
        key2idx = {}
        def make_key(pdb_id, chain_id, ref_residue, position, mutated_residue):
            return (pdb_id, chain_id, ref_residue, position, mutated_residue)

        for idx, sample in self.data.items():
            key = make_key(sample["pdb_id"], sample["chain_id"], sample["ref_residue"], sample["resseq_position"], sample["mutated_residue"])
            assert(key not in key2idx)
            key2idx[key] = idx

        from pandas import read_csv
        csv_data = read_csv(csv_path, sep='\t')
        for _, sample in csv_data.iterrows():
            pdb_id = sample["PDB Id"].lower()
            chain_id = sample["Mutated Chain"]
            if pdb_id == "1rtp":
                chain_id = '1'

            mutation = sample["Mutation_PDB"]
            ref_residue = mutation[0]
            position = int(mutation[1:-1])
            mutated_residue = mutation[-1]

            is_in_309 = True if sample["Label"] == 1 else False
            uniprot_id = sample["UniProt"]
            if pdb_id in ["1yyj", "2imm", "2a01"]:
                uniprot_id = '-'

            key = make_key(pdb_id, chain_id, ref_residue, position, mutated_residue)
            idx = key2idx[key]

            self.data[idx]["is_in_309"] = is_in_309
            self.data[idx]["uniprot_id"] = uniprot_id

    def export(self, csv_path: str) -> None:
        super().export(csv_path, columns=self.S350_COLUMNS_LIST)
