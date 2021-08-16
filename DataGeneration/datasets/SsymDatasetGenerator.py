from .DatasetGenerator import DatasetGenerator

class SsymDatasetGenerator(DatasetGenerator):
    SSYM_COLUMNS_LIST=["id", "pdb_id", "chain_id", "uniprot_id", "ref_residue", "seq_position", "resseq_position", "mutated_residue", "pH", "T", "ddG", "sequence"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def __make_key(idx, label: str) -> None:
        label = "F" if label == "forward" else "R"
        key = str(idx) + label
        return key

    def read_Ssym(self, csv_path: str) -> None:
        from pandas import read_csv
        csv_data = read_csv(csv_path)
        csv_data.drop(index=0, inplace=True) # remove extra header
        for idx, sample in csv_data.iterrows():
            assert(sample["Wild Type"].upper() == sample["Mutant.1"].upper())
            assert(sample["Wild Type.1"].upper() == sample["Mutant"].upper())

            def add_sample(label, pdb_id, chain_id, ref_residue, resseq_position, mutated_residue, pH, T, ddG):               
                from Bio.SeqUtils import seq1
                ref_residue = seq1(ref_residue)
                mutated_residue = seq1(mutated_residue)            

                key = self.__make_key(idx, label)
                self.add(key,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    ref_residue=ref_residue,
                    resseq_position=resseq_position,
                    mutated_residue=mutated_residue,
                    pH=pH,
                    T=T,
                    ddG=ddG
                )

            pH = float(sample["Condition pH"].replace(',', '.'))
            T = float(sample["Condition T"].replace(',', '.'))
            ddGExp = float(sample["ΔΔGexp"].replace(',', '.'))

            assert(sample["Res Num"] == sample["Res Number"]) # optional but good to know that this is the case
            add_sample("forward", sample["PDB"], sample["Chain"], sample["Wild Type"], sample["Res Num"], sample["Mutant"],  pH, T, ddGExp)
            add_sample("reverse", sample["PDB.1"], sample["Chain.1"], sample["Wild Type.1"], sample["Res Number"], sample["Mutant.1"], pH, T, -ddGExp)

    def read_Ssym_PremPS(self, csv_path: str) -> None:
        pdb2uniprot = {}
        import pandas as pd
        df = pd.read_csv(csv_path, sep='\t')
        df["PDB Id"] = df["PDB Id"].str.lower()
        df["Extended PDB Id"] = df["PDB Id"] + df["Mutated Chain"]
        df = df[["Extended PDB Id", "UniProt"]]
        series = pd.Series(df["UniProt"].values, index=df["Extended PDB Id"].values)
        pdb2uniprot = series.to_dict()
        for key, sample in self.data.items():
            extended_pdb_id = sample["pdb_id"] + sample["chain_id"]
            self.data[key]["uniprot_id"] = pdb2uniprot[extended_pdb_id]

    def validate_data(self, resseq_dump_dir_path):
        super().validate_data(resseq_dump_dir_path)
        forward_samples = { key[:-1] : sample for key, sample in self.data.items() if sample["id"][-1] == "F"}
        reverse_samples = { key[:-1] : sample for key, sample in self.data.items() if sample["id"][-1] == "R"}
        for key, fwd_sample in forward_samples.items():
            rev_sample = reverse_samples[key]
            assert(fwd_sample["ref_residue"] == rev_sample["mutated_residue"])
            assert(fwd_sample["mutated_residue"] == rev_sample["ref_residue"])
            assert(fwd_sample["ddG"] == -rev_sample["ddG"])
            assert(fwd_sample["pH"] == rev_sample["pH"])
            assert(fwd_sample["T"] == rev_sample["T"])

    def export(self, samples_csv_path, idx_csv_path):
        super().export(samples_csv_path, columns=self.SSYM_COLUMNS_LIST)

        forward_samples = { key[:-1] : sample for key, sample in self.data.items() if sample["id"][-1] == "F"}
        reverse_samples = { key[:-1] : sample for key, sample in self.data.items() if sample["id"][-1] == "R"}
        assert(forward_samples.keys() == reverse_samples.keys())
        with open(idx_csv_path, "w") as f:
            for key in forward_samples.keys():
                f.write("%s,%s\n" % (key + 'F', key + 'R'))

