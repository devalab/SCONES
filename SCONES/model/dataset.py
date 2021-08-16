import torch
import numpy as np

from torch.utils.data import Dataset
from utils.constants import AA_ID_DICT, AA_SASA

class PreprocessedDataset(Dataset):
    def __init__(self, radius_cutoff, device, logger=None):
        self.logger = logger
        if logger is None:
            import logging
            self.logger = logging
        self.device = device
        self.reset(radius_cutoff)

    def reset(self, radius_cutoff):
        self.dataset = []
        self.structure_cache = {}
        self.radius_cutoff = radius_cutoff
        
        self.counts = {
            "no_target_skipped" : 0,
            "insufficient_structure" : 0,
            "insufficient_env_structure" : 0,
            "no_struct_at_mutation" : 0,
            "ddG_out_of_range" : 0,
            "proline_mutation" : 0,
            "not_enough_neighbors" : 0,
            "missing_nonstd_residue" : 0,
            "pH_out_of_range_skipped" : 0
        }

    def add_dataset(self, dataset, training_filters):
        for sample in dataset:
            self.add(sample, training_filters)

    def preprocess_structure(self, structure):
        primary = list(structure["sequence"])
        seqlen = len(primary)

        # [seqlen]
        transformed_primary = np.vectorize(AA_ID_DICT.get)(primary)
        transformed_primary = transformed_primary
        transformed_primary = torch.from_numpy(transformed_primary)

        assert(structure["distance_map_cb"].shape[0] == seqlen)
        
        mask = structure["structure_mask"]
        pos_avail = np.count_nonzero(mask)
        if pos_avail / len(mask) < 0.5:
            self.counts["insufficient_structure"] += 1
            return None

        sas = torch.from_numpy(structure["sas"])
        sas = sas.float()
            
        processed_struct = {
            "primary" : primary,
            "primary_idx" : transformed_primary,
            "sas" : sas,
            "mask" : mask,
            "dmap_cb" : structure["distance_map_cb"],
            "omap_w" : structure["orientation_map_omega"],
            "omap_theta" : structure["orientation_map_theta"],
            "omap_phi" : structure["orientation_map_phi"],
        }
        return processed_struct

    def add(self, sample, training_filters=False):
        def add_sample(id, dataset_sample, sequence, ref_residue, mutated_residue, primary_rt, primary_props_rt, primary_mt, primary_props_mt, position, neighbors, edge_features, edge_features_rev, target):
            sample = {
                "id" : id,
                "dataset_sample" : dataset_sample,
                "sequence" : sequence,
                "ref_residue" : ref_residue,
                "mutated_residue" : mutated_residue,
                "primary_rt" : primary_rt,
                "primary_props_rt" : primary_props_rt,
                "primary_mt" : primary_mt,
                "primary_props_mt" : primary_props_mt,
                "position"  : position,
                "neighbors" : neighbors,
                "edge_features" : edge_features,
                "edge_features_rev" : edge_features_rev,
                "target" : target,
            }

            if self.device:
                for key, value in sample.items():
                    if torch.is_tensor(value):
                        sample[key] = sample[key].to(self.device)

            self.dataset.append(sample)                

        target_ddG = sample["ddG"]
        if np.isnan(target_ddG):
            self.counts["no_target_skipped"] += 1
            return

        if training_filters and np.abs(target_ddG) > 8.0:
            self.counts["ddG_out_of_range"] += 1
            return

        pH = sample["dataset_sample"].get("pH", 7.0)
        if training_filters and (2.0 > pH or pH > 12.0):
            self.counts["pH_out_of_range_skipped"] += 1
            return

        if 'X' in sample["structure"]["sequence"]:
            self.counts["missing_nonstd_residue"] += 1
            return

        raw_structure = sample["structure"]
        idx = raw_structure["id"]
        if idx not in self.structure_cache:
            retval = self.preprocess_structure(raw_structure)
            if retval is None:
                return
            self.structure_cache[idx] = retval
        processed_structure = self.structure_cache[idx]

        sequence = processed_structure["primary"]
        primary = processed_structure["primary_idx"]
        mask = processed_structure["mask"]
        position_c = sample["position"] # central residue position

        sas = processed_structure["sas"]
        dmap_cb = processed_structure["dmap_cb"]
        omap_w = processed_structure["omap_w"]
        omap_theta = processed_structure["omap_theta"]
        omap_phi = processed_structure["omap_phi"]

        # def is_struct_present(position):
        #     if np.isnan(sas[position]):
        #         return False
        #     if np.any(np.isnan(dmap_cb[position, :])):
        #         if mask[position] & 5 != 5:
        #             return False
        #     elif mask[position] & 7 != 7:
        #         return False
        #     return True

        def is_struct_present(position, sas=True):
            if sas and mask[position] & 16 == 0:
                return False
            if primary[position] == AA_ID_DICT['G']:
                if mask[position] & 5 != 5:
                    return False
            elif mask[position] & 7 != 7:
                return False
            return True

        if not is_struct_present(position_c):
            self.counts["no_struct_at_mutation"] += 1
            return

        # find neighbors using Cb
        cb_dist = dmap_cb[position_c, :]
        neighbors = np.where((cb_dist > 0) & (cb_dist < self.radius_cutoff))[0]
        num_neighbors = len(neighbors)
        assert (num_neighbors > 0)

        if training_filters and num_neighbors < 5:
            self.counts["not_enough_neighbors"] += 1
            return

        if training_filters and not np.all([is_struct_present(i) for i in neighbors]):
            self.counts["insufficient_env_structure"] += 1
            return

        neighbors = neighbors[[is_struct_present(i) for i in neighbors]]
        num_neighbors = len(neighbors)
        
        edge_features = torch.zeros((num_neighbors, 16))
        for i in range(num_neighbors):
            position_r = neighbors[i]

            cb_dist = dmap_cb[position_c, position_r]
            w_ang = omap_w[position_c, position_r]
            
            theta_ang12, theta_ang21 = omap_theta[position_c, position_r], omap_theta[position_r, position_c]
            phi_ang12, phi_ang21 = omap_phi[position_c, position_r], omap_phi[position_r, position_c]

            means8 = np.asarray([6.2193589e+00, 1.6547333e-01, 2.8232558e-02, 4.9734241e-03, 2.9435262e-01, 1.3330153e-01])
            scales8 = np.asarray([1.0191348, 0.02915691, 0.01040565, 0.00290022, 0.06066265, 0.0461563])

            edge_features[i, 0] = cb_dist
            edge_features[i, 1] = 1/cb_dist
            edge_features[i, 2] = (1/cb_dist)**2
            edge_features[i, 3] = (1/cb_dist)**3
            edge_features[i, 4] = np.exp(-cb_dist/5)
            edge_features[i, 5] = np.exp(-cb_dist/3)
            edge_features[i, :6] -= means8
            edge_features[i, :6] /= scales8

            if not np.isnan(w_ang):
                edge_features[i, 6] = np.sin(w_ang)
                edge_features[i, 7] = np.cos(w_ang)

            if not np.isnan(theta_ang12):
                edge_features[i, 8] = np.sin(theta_ang12)
                edge_features[i, 9] = np.cos(theta_ang12)
            
            if not np.isnan(theta_ang21):
                edge_features[i, 10] = np.sin(theta_ang21)
                edge_features[i, 11] = np.cos(theta_ang21)

            if not np.isnan(phi_ang12):
                edge_features[i, 12] = np.sin(phi_ang12)
                edge_features[i, 13] = np.cos(phi_ang12)

            if not np.isnan(phi_ang21):
                edge_features[i, 14] = np.sin(phi_ang21)
                edge_features[i, 15] = np.cos(phi_ang21)

        edge_features_rev = edge_features.clone()
        for i in range(num_neighbors):
            edge_features_rev[i, 8] = edge_features[i, 10]
            edge_features_rev[i, 9] = edge_features[i, 11]
            edge_features_rev[i, 10] = edge_features[i, 8]
            edge_features_rev[i, 11] = edge_features[i, 9]
            edge_features_rev[i, 12] = edge_features[i, 14]
            edge_features_rev[i, 13] = edge_features[i, 15]
            edge_features_rev[i, 14] = edge_features[i, 12]
            edge_features_rev[i, 15] = edge_features[i, 13]

        ref_residue = primary[position_c].item()
        mutated_residue = AA_ID_DICT[sample["mutated_residue"]]

        primary_rt = primary.clone()
        primary_rt[position_c] = ref_residue

        primary_mt = primary.clone()
        primary_mt[position_c] = mutated_residue

        target = torch.zeros(1)
        target[0] = target_ddG

        primary_props_rt = sas.clone().reshape(1, -1)
        primary_props_mt = sas.clone().reshape(1, -1)
        primary_props_mt[0, position_c] *= AA_SASA[sample["mutated_residue"]]/AA_SASA[sample["ref_residue"]] # heuristic approximation of mutated residue SAS
        add_sample(sample["id"], sample["dataset_sample"], sequence, sample["ref_residue"], sample["mutated_residue"], primary_rt, primary_props_rt, primary_mt, primary_props_mt, position_c, neighbors, edge_features, edge_features_rev, target)

    def summary(self):
        self.logger.info("Dropped %d samples as they did not have ddG" % self.counts["no_target_skipped"])
        self.logger.info("Dropped %d samples as they have insufficient structure" % self.counts["insufficient_structure"])
        self.logger.info("Dropped %d samples as they have missing neighbor structure" % self.counts["insufficient_env_structure"])
        self.logger.info("Dropped %d samples as their pH during measurement was not in acceptable range" % self.counts["pH_out_of_range_skipped"])
        self.logger.info("Dropped %d samples as the mutation site didn't have local structure information" % self.counts["no_struct_at_mutation"])
        self.logger.info("Dropped %d samples as target was out of range" % self.counts["ddG_out_of_range"])
        self.logger.info("Dropped %d samples as mutation involved proline" % self.counts["proline_mutation"])
        self.logger.info("Dropped %d samples as they had too few neighbors" % self.counts["not_enough_neighbors"])
        self.logger.info("Dropped %d samples as they had missing or non-standard residue" % self.counts["missing_nonstd_residue"])

        # print("Average number of neighbors:", np.mean(neighbor_counts))
        # print("Average edge features:", np.mean(edge_features_pool, axis=0), np.std(edge_features_pool, axis=0))
        stabalizing_samples, neutral_samples, destabalizing_samples = 0, 0, 0
        for sample in self.dataset:
            target_ddG = sample["target"].item()
            if target_ddG > 0.5: stabalizing_samples += 1
            elif target_ddG < -0.5: destabalizing_samples += 1
            else: neutral_samples += 1
        self.logger.info("Stabalizing samples: %d, neutral samples: %d, destabalizing samples: %d" % (stabalizing_samples, neutral_samples, destabalizing_samples))

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        sample = self.dataset[idx]
        return sample

    
