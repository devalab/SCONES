import typing
import numpy as np

from utils.sample_key import SampleKey

class MergeDuplicates:
    def __init__(self):
        from collections import defaultdict
        self.data = {}
        self.sample_groups = defaultdict(list)

    def merge(self, data: dict, resseq_dump_dir_path: str, default_pH: float=7.0, default_T: float=0):
        self.data = data
        for idx, sample in data.items():
            pdb_id = sample["pdb_id"]
            chain_id = sample["chain_id"]
            extended_pdb_id = pdb_id + chain_id
            ref_residue = sample["ref_residue"]
            mutated_residue = sample["mutated_residue"]
            if "sequence" not in sample or "seq_position" not in sample:
                import os
                fasta_path = os.path.join(resseq_dump_dir_path, extended_pdb_id +  ".fasta")
                resseq_path = os.path.join(resseq_dump_dir_path, extended_pdb_id + ".resseq")

                from Bio import SeqIO
                sequence = next(SeqIO.parse(fasta_path, "fasta"))
                sequence = str(sequence.seq)

                resseq = next(SeqIO.parse(resseq_path, "fasta"))
                resseq = str(resseq.seq)
                resseq = [int(i) for i in resseq.split(',')]
                assert(len(resseq) == len(sequence))

                resseq = {k : v for v, k in enumerate(resseq)}

                resseq_position = sample["resseq_position"]
                assert(resseq_position in resseq)
                seq_position = resseq[resseq_position]
            else:
                sequence = sample["sequence"]
                seq_position = sample["seq_position"]

            assert(sequence[seq_position] == ref_residue)

            pH = sample["pH"] if "pH" in sample else default_pH
            T = sample["T"] if "T" in sample else default_T

            key = SampleKey(sequence, ref_residue, seq_position, mutated_residue, pH, T)
            self.sample_groups[key].append(idx)

    def export(self, csv_path: str, mergefn, columns: typing.Optional[list]=None) -> None:
        records = {}        
        for group in self.sample_groups.values():
            samples = [self.data[idx] for idx in group]
            merged_sample = mergefn(samples)
            if merged_sample is not None:
                id = merged_sample["id"]
                records[id] = merged_sample
        
        import pandas as pd
        df = pd.DataFrame.from_dict(records).T
        if columns:
            for column in columns:
                assert(column in df.columns)
            df = df.reindex(columns, axis=1)
            assert(df.columns.size == len(columns))
        df.to_csv(csv_path, index=False)

