import torch
import os
import numpy as np
import pickle

class SCONESDataset:
    def __init__(self):
        self.mutation_data = []
        self.chain_data = {}

        self.pdb_cache = {}
        self.dssp_cache = {}
        self.discarded_pdb_chain_ids = []

        self.opt = {
            "max_length" : 3000,

            "aligner_config" : {
                "identity_score" : 5,
                "mismatch_penalty" : -4,
                "gap_init_penalty" : -5,
                "gap_extend_penalty" : -1
            },
        }

    def summary(self):
        print("total mutations: ", len(self.mutation_data))
        print("total polypeptides: ", len(self.chain_data))
        print("data memory size: %d KB" % (len(pickle.dumps(self.mutation_data)) / 1024))
        print("polypeptides memory size: %d KB" % (len(pickle.dumps(self.chain_data)) / 1024))

    def load_dataset(self, csv_path, PDB_PATH):
        assert(os.path.isdir(PDB_PATH))

        from pandas import read_csv
        dataset = read_csv(csv_path)
        for _, sample in dataset.iterrows():
            status, errmsg = self.add_sample(sample, PDB_PATH)
            if status is False:
                print("[SKIP]", sample["id"], " >> errmsg:", errmsg, flush=True)

    def add_sample(self, sample, PDB_PATH):
        print("Processing sample with id", sample["id"])
        pdb_id = sample["pdb_id"]
        # if not isinstance(pdb_id, str):
        #     return False, "no pdb"

        ref_sequence = sample["sequence"]
        if len(ref_sequence) > self.opt["max_length"]:
            return False, "polypeptide exceeds maximum size limits (ref: %d)" % len(ref_sequence)

        chain_id = sample["chain_id"]
        extended_pdb_id = pdb_id + chain_id

        if extended_pdb_id in self.discarded_pdb_chain_ids:
            return False, "polypeptide discared previously"

        if pdb_id not in self.pdb_cache:
            pdb_path = os.path.join(PDB_PATH, pdb_id + ".cif")
            print("[CIF PARSE] parsing", pdb_path, flush=True)
            from Bio.PDB.MMCIFParser import FastMMCIFParser
            mmcif_parser = FastMMCIFParser(QUIET=True)
            struct = mmcif_parser.get_structure(pdb_id, pdb_path)
            print("[DSSP PARSE] parsing", pdb_path, flush=True)
            try:
                from Bio.PDB.DSSP import dssp_dict_from_pdb_file
                dssp_info = dssp_dict_from_pdb_file(pdb_path)
            except:
                return False, "DSSP failed"
            self.pdb_cache[pdb_id] = struct
            self.dssp_cache[pdb_id] = dssp_info

        if extended_pdb_id not in self.chain_data:
            print("[POLY PARSE] parsing", extended_pdb_id, flush=True)
            polypeptide = self.__parse_polypeptide(self.pdb_cache[pdb_id], self.dssp_cache[pdb_id], chain_id)
            if polypeptide is None:
                print("[SKIP] discarded", extended_pdb_id, flush=True)
                self.discarded_pdb_chain_ids.append(extended_pdb_id)
                return False, "failed to process structure"
            polypeptide["id"] = extended_pdb_id
            self.chain_data[extended_pdb_id] = polypeptide

        polypeptide = self.chain_data[extended_pdb_id]
        sequence = polypeptide["sequence"]
        position = sample["seq_position"]

        if sequence != ref_sequence:
            from Bio import pairwise2
            aligner_config = self.opt["aligner_config"]
            alignments = pairwise2.align.globalms(
                ref_sequence, sequence,
                aligner_config["identity_score"], aligner_config["mismatch_penalty"], aligner_config["gap_init_penalty"], aligner_config["gap_extend_penalty"],
                penalize_end_gaps=False,
                one_alignment_only=True
            )

            alignment = alignments[0]
            seqA = alignment.seqA
            seqB = alignment.seqB
            seqlen = len(seqA)
            assert(len(seqA) == len(seqB))
            assert(alignment.start == 0 and alignment.end == seqlen)

            mismatch = False
            posA, posB = 0, 0
            for i in range(len(seqA)):
                if seqA[i] != '-':
                    posA += 1
                if seqB[i] != '-':
                    posB += 1
                if posA == position + 1:
                    mismatch = seqA[i] != seqB[i]
                    break
            
            if posB == 0 or mismatch:
                return False, "no structure at the mutation site or a mismatch"

            print(seqA)
            print(seqB)
            print(position + 1, posA, posB, sequence[posB - 1], sample["ref_residue"])
            position = posB - 1

        assert(sequence[position] == sample["ref_residue"])

        if sequence[position] != sample["ref_residue"]:
            print("residue mismatch >> found: %c, expected: %c, position: %d" % (sequence[position], sample["ref_residue"], position))
            print(sequence)
            return False, "residue mismatch"
        
        self.mutation_data.append({
            "id" : sample["id"],
            "structure" : polypeptide,
            "ref_residue" : sample["ref_residue"],
            "position" : position,
            "mutated_residue" : sample["mutated_residue"],
            "ddG" : sample["ddG"],
            "dataset_sample" : sample
        })

        return True, None

    def __parse_polypeptide(self, struct, dssp_info, chain_id):
        model = struct[0]

        from Bio.SeqUtils import seq1
        try:
            chain = model[chain_id]
        except:
            print("[SKIP] chain %s does not exist in model for" % chain_id, struct.get_id())
            return None

        from Bio.PDB.Polypeptide import is_aa
        residues = [r for r in chain.get_residues() if is_aa(r)]
        resseq_ids = [r.get_id()[1] for r in residues]
        sequence = [seq1(r.get_resname()) for r in residues]
        sequence = "".join(sequence)
        assert(len(residues) == len(resseq_ids))
        assert(len(residues) == len(sequence))
        num_residues = len(residues)

        structure_mask = np.zeros(num_residues, dtype=np.int)
        CA_MASK, CB_MASK, N_MASK, SS_MASK, SAS_MASK = 1, 2, 4, 8, 16
        CA_coords = [None] * num_residues
        CB_coords = [None] * num_residues
        N_coords = [None] * num_residues
        sstruct = np.full((num_residues), '-', dtype=np.str)
        sas = np.full((num_residues), np.nan, dtype=np.float)

        def is_atom_present(mask, keys):
            return mask & keys == keys

        for i, residue_i in enumerate(residues):
            if "CA" in residue_i:
                structure_mask[i] |= CA_MASK
                CA_coords[i] = residue_i["CA"].coord
            if "CB" in residue_i:
                structure_mask[i] |= CB_MASK
                CB_coords[i] = residue_i["CB"].coord
            if "N" in residue_i:
                structure_mask[i] |= N_MASK
                N_coords[i] = residue_i["N"].coord

        dmap_cb = np.full((num_residues, num_residues), np.nan, dtype=np.float)
        omap_omega = np.full((num_residues, num_residues), np.nan, dtype=np.float)
        omap_theta = np.full((num_residues, num_residues), np.nan, dtype=np.float)
        omap_phi = np.full((num_residues, num_residues), np.nan, dtype=np.float)
        sstruct = np.full((num_residues), '-', dtype=np.str)
        sas = np.full((num_residues), np.nan, dtype=np.float)

        dssp, dssp_keys = dssp_info
        for i, residue_i in enumerate(residues):
            dssp_idx = (chain_id, residue_i.get_id())
            if dssp_idx in dssp_keys:
                dssp_i = dssp[dssp_idx]
                sstruct[i] = dssp_i[1]
                structure_mask[i] |= SS_MASK 
                if dssp_i[2] != 'NA':
                    structure_mask[i] |= SAS_MASK
                    sas[i] = float(dssp_i[2])
                dssp_aa = dssp_i[0]
                if str.islower(dssp_aa):
                    dssp_aa = 'C'
                assert(dssp_aa == sequence[i])

        for i, residue_i in enumerate(residues): 
            for j, residue_j in enumerate(residues):
                from Bio.PDB.vectors import Vector, calc_dihedral, calc_angle
                if i == j:
                    dmap_cb[i][j] = np.nan
                    omap_omega[i][j] = np.nan
                    omap_theta[i][j] = np.nan
                    omap_phi[i][j] = np.nan
                    continue

                Cb_i = CB_coords[i] if is_atom_present(structure_mask[i], CB_MASK) else None
                if sequence[i] == 'G':
                    Cb_i = CA_coords[i] if is_atom_present(structure_mask[i], CA_MASK) else None

                Cb_j = CB_coords[j] if is_atom_present(structure_mask[j], CB_MASK) else None
                if sequence[j] == 'G':
                    Cb_j = CA_coords[j] if is_atom_present(structure_mask[j], CA_MASK) else None
    
                # symmetric features
                if i > j:
                    if Cb_i is not None and Cb_j is not None:
                        diff = Cb_i - Cb_j
                        dmap_cb[i][j] = dmap_cb[j][i] = np.linalg.norm(diff)

                    if is_atom_present(structure_mask[i], CA_MASK | CB_MASK) and is_atom_present(structure_mask[j], CA_MASK | CB_MASK):
                        v1 = Vector(CA_coords[i])
                        v2 = Vector(CB_coords[i])
                        v3 = Vector(CB_coords[j])
                        v4 = Vector(CA_coords[j])
                        omap_omega[i][j] = omap_omega[j][i] = calc_dihedral(v1, v2, v3, v4)

                # asymmetric features
                if is_atom_present(structure_mask[i], CA_MASK | CB_MASK | N_MASK) and Cb_j is not None:
                    v1 = Vector(N_coords[i])
                    v2 = Vector(CA_coords[i])
                    v3 = Vector(CB_coords[i])
                    v4 = Vector(Cb_j)
                    omap_theta[i][j] = calc_dihedral(v1, v2, v3, v4)

                if is_atom_present(structure_mask[i], CA_MASK | CB_MASK) and Cb_j is not None:
                    v1 = Vector(CA_coords[i])
                    v2 = Vector(CB_coords[i])
                    v3 = Vector(Cb_j)
                    omap_phi[i][j] = calc_angle(v1, v2, v3)

        for i in range(num_residues):
            for j in range(num_residues):
                if i == j:
                    continue
                if is_atom_present(structure_mask[i], CB_MASK) and is_atom_present(structure_mask[j], CB_MASK):
                    assert(not np.isnan(dmap_cb[i][j]))
                if is_atom_present(structure_mask[i], CA_MASK | CB_MASK) and is_atom_present(structure_mask[j], CA_MASK | CB_MASK) & 3 == 3:
                    assert(not np.isnan(omap_omega[i][j]))
                if is_atom_present(structure_mask[i], CA_MASK | CB_MASK | N_MASK) and is_atom_present(structure_mask[j], CB_MASK):
                    assert(not np.isnan(omap_theta[i][j]))
                if is_atom_present(structure_mask[i], CA_MASK | CB_MASK) and is_atom_present(structure_mask[j], CB_MASK):
                    assert(not np.isnan(omap_phi[i][j]))

        return {
            "id" : struct.get_id(),
            "sequence" : sequence,
            "resseq_ids" : resseq_ids,
            "structure_mask" : structure_mask,
            "secondary_structure" : sstruct,
            "sas" : sas,
            "distance_map_cb" : dmap_cb,
            "orientation_map_omega" : omap_omega,
            "orientation_map_theta" : omap_theta,
            "orientation_map_phi" : omap_phi
        }

    def serialize(self, filename):
        with open(filename, "wb") as file:
            pickle.dump((self.chain_data, self.mutation_data), file)

    def deserialize(self, filename):
        with open(filename, "rb") as file:
            self.chain_data, self.mutation_data = pickle.load(file)   

    def __len__(self):
        return len(self.mutation_data)
    
    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        sample = self.mutation_data[idx]
        return sample