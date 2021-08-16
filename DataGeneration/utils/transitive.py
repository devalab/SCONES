import typing
import numpy as np

from collections import defaultdict
from utils.sample_key import SampleKey, ALIGNMENT_SCORE_MIN
from utils.get_aligned_position import get_aligned_position

class TransitiveClosure:
    def __init__(self):
        self.source_data = {}
        self.source_key2id = {}
        self.data = {}
        self.transitive_pairs = set()

    def load_dataset(self, data: dict, resseq_dump_dir_path: str, default_pH: float=7.0, default_T: float=0):
        self.source_data = data
        for idx, sample in enumerate(data.values()):
            pdb_id, chain_id = sample["pdb_id"], sample["chain_id"]
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

            pH = sample["pH"] if "pH" in sample and not np.isnan(sample["pH"]) else default_pH
            T = sample["T"] if "T" in sample and not np.isnan(sample["T"]) else default_T

            new_sample = {
                "sequence" : sequence,
                "seq_position" : seq_position,
                "ref_residue" : ref_residue,
                "mutated_residue" : mutated_residue,
                "pH" : pH,
                "T" : T,
            }

            key = self.get_key(new_sample)
            if key in self.data:
                continue
            self.data[key] = new_sample
            self.source_key2id[key] = sample["id"]

    def get_key(self, sample):
        primary_rp = sample["sequence"]
        ref_residue = sample["ref_residue"]
        seq_position = sample["seq_position"]
        mutated_residue = sample["mutated_residue"]
        pH = sample["pH"]
        T = sample["T"]
        assert(primary_rp[seq_position] == ref_residue)
        return SampleKey(primary_rp, ref_residue, seq_position, mutated_residue, pH, T)

    def get_transitive_key(self, key1: SampleKey, key2: SampleKey, reflexive_result_allowed: bool=False):
        primary_rp1 = key1.primary_rp
        ref_residue1, position1, mutated_residue1 = key1.mutation
        pH1, T1 = key1.pH, key1.T

        primary_rp2 = key2.primary_rp
        ref_residue2, position2, mutated_residue2 = key2.mutation
        pH2, T2 = key2.pH, key2.T

        if not (mutated_residue1 == ref_residue2):
            return None

        if not SampleKey.is_pH_same(pH1, pH2):
            return None

        if not SampleKey.is_T_same(T1, T2):
            return None

        if not reflexive_result_allowed:
            if ref_residue1 == mutated_residue2: # reflexive key
                return None

        primary_mp1 = list(primary_rp1)
        primary_mp1[position1] = mutated_residue1
        primary_mp1 = "".join(primary_mp1)

        if primary_mp1 == primary_rp2:
            if position1 != position2:
                return None
        else:
            score, new_position2 = get_aligned_position(primary_mp1, position1, primary_rp2)
            if score is None or score < ALIGNMENT_SCORE_MIN:
                return None
            assert(primary_mp1[position1] == primary_rp2[position2])
            if new_position2 != position2:
                return None
        
        return SampleKey(primary_rp1, ref_residue1, position1, mutated_residue2, (pH1 + pH2) * 0.5, (T1 + T2) * 0.5)

    def transitive_closure(self, max_iters: int=1):
        def make_transitive_sample(sample1, sample2):
            assert(sample1["mutated_residue"] == sample2["ref_residue"])
            seq_position = sample1["seq_position"]
            sample = {
                "sequence" : sample1["sequence"],
                "seq_position" : seq_position,
                "ref_residue" : sample1["ref_residue"],
                "mutated_residue" : sample2["mutated_residue"],
                "pH" : (sample1["pH"] + sample2["pH"])/2,
                "T" : (sample1["T"] + sample2["T"])/2,
            }
            assert(sample["sequence"][seq_position] == sample["ref_residue"])
            return sample

        def assert_transitivity(sample1, sample2, sample3):
            assert(sample1["sequence"][sample1["seq_position"]] == sample1["ref_residue"])
            assert(sample2["sequence"][sample2["seq_position"]] == sample2["ref_residue"])
            assert(sample3["sequence"][sample3["seq_position"]] == sample3["ref_residue"])
            assert(sample1["ref_residue"] == sample3["ref_residue"])
            assert(sample1["mutated_residue"] == sample2["ref_residue"])
            assert(sample2["mutated_residue"] == sample3["mutated_residue"])

        def assert_key(key, sample):
            assert(self.get_key(sample) == key)

        for itr in range(max_iters):
            current_size = len(self.data)
            print("[iter %d] current size:" % itr, current_size)

            new_samples = {}
            for idx1, (key1, sample1) in enumerate(self.data.items()):
                for key2, sample2 in self.data.items():
                    transitive_key = self.get_transitive_key(key1, key2, False) # without reflexive result
                    if transitive_key and transitive_key not in self.data:
                        self.transitive_pairs.add((key1, key2, transitive_key))
                        sample3 = make_transitive_sample(sample1, sample2)
                        new_samples[transitive_key] = sample3

            for key, sample in new_samples.items():
                self.data[key] = sample

            for key1, key2, key3 in self.transitive_pairs:
                sample1 = self.data[key1]
                sample2 = self.data[key2]
                sample3 = self.data[key3]
                assert_key(key1, self.data[key1])
                assert_key(key2, self.data[key2])
                assert_key(key3, self.data[key3])
                assert_transitivity(sample1, sample2, sample3)

            new_size = len(self.data)
            if current_size == new_size:
                break

        print("transitive pairs:", len(self.transitive_pairs))

        for key1, key2, key3 in self.transitive_pairs:
            sample1 = self.data[key1]
            sample2 = self.data[key2]
            sample3 = self.data[key3]
            #print(sample1["ref_residue"], "->", sample1["mutated_residue"], sample1["pH"], ";", sample2["ref_residue"], "->", sample2["mutated_residue"], sample2["pH"], ";", sample3["ref_residue"], "->", sample3["mutated_residue"], sample3["pH"])
            assert(sample1["ref_residue"] == sample3["ref_residue"])
            assert(sample1["mutated_residue"] == sample2["ref_residue"])
            assert(sample2["mutated_residue"] == sample3["mutated_residue"])
    
    def export(self, samples_csv_path: str, idx_csv_path: str, mergefn, columns: typing.Optional[list]=None):
        import copy
        records = copy.deepcopy(self.source_data)
        key2id = copy.deepcopy(self.source_key2id)
        indices = []        
        for key1, key2, key3 in self.transitive_pairs:
            id1, id2 = key2id[key1], key2id[key2]
            sample1, sample2 = records[id1], records[id2]
            sample3 = mergefn(sample1, sample2)
            id3 = sample3["id"]
            key2id[key3] = id3
            records[id3] = sample3
            indices.append((id1, id2, id3))

        import pandas as pd
        df = pd.DataFrame.from_dict(records).T
        if columns:
            for column in columns:
                assert(column in df.columns)
            df = df.reindex(columns, axis=1)
            assert(df.columns.size == len(columns))
        df.to_csv(samples_csv_path, index=False)
        
        with open(idx_csv_path, "w") as f:
            for key1, key2, key3 in self.transitive_pairs:
                id1, id2, id3 = key2id[key1], key2id[key2], key2id[key3]
                f.write("%s,%s,%s\n" % (id1, id2, id3))

