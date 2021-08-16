import torch
import torch.nn as nn
import numpy as np

from .base import BaseNetwork

class EmbeddingNetwork(nn.Module):
    def __init__(self, out_features):
        super(EmbeddingNetwork, self).__init__()
        self.primary_embedding = nn.Embedding(20, out_features)
        self.init_embeddding()

    def init_embeddding(self):
        from quantiprot.utils.feature import Feature
        from quantiprot.utils.mapping import Mapping
        from quantiprot.metrics.aaindex import get_aa2hydropathy, get_aa2volume, get_aaindex_file

        def normalize(feature, mean=None, scale=None):
            mapping = feature.function.mapping
            if mean is None or scale is None:
                values = list(mapping.values())
                mean = np.mean(values)
                scale = np.std(values)
            mapping.update({key: (mapping[key] - mean) / scale for key in mapping.keys()})
            return Feature(Mapping(feature.name, mapping))

        aa2charge = Mapping("formal_charge", {'A': 0, 'C': 0, 'D': -1, 'E': -1,
                                              'F': 0, 'G': 0, 'H': 1, 'I': 0,
                                              'K': 1, 'L': 0, 'M': 0, 'N': 0,
                                              'P': 0, 'Q': 0, 'R': 1, 'S': 0,
                                              'T': 0, 'V': 0, 'W': 0, 'Y': 0,
                                             }, default=None)
    
        f_steric_param = Feature(get_aaindex_file("CHAM810101")) # steric parameter
        f_polarity = Feature(get_aaindex_file("GRAR740102")) # polarity 
        f_ras_tripep = Feature(get_aaindex_file("CHOC760101")) # residue accessible surface area in tripeptide

        from utils.constants import AA_ID_DICT
        embedding_weights = torch.zeros((20, 6), dtype=torch.float32)
        reverse_id = { v : k for k, v in AA_ID_DICT.items() }
        for idx, aa in reverse_id.items():
            embedding_weights[idx, 0] = aa2charge.mapping[aa]
            embedding_weights[idx, 1] = normalize(Feature(get_aa2volume())).function.mapping[aa]
            embedding_weights[idx, 2] = normalize(Feature(get_aa2hydropathy())).function.mapping[aa]
            embedding_weights[idx, 3] = normalize(f_steric_param).function.mapping[aa]
            embedding_weights[idx, 4] = normalize(f_polarity).function.mapping[aa]
            embedding_weights[idx, 5] = normalize(f_ras_tripep).function.mapping[aa]

        self.primary_embedding.weight = torch.nn.Parameter(embedding_weights)

    def disable_training(self):
        self.primary_embedding.weight.requires_grad = False

    def enable_training(self):
        self.primary_embedding.weight.requires_grad = True

    def forward(self, primary):
        primary = primary.reshape(1, -1)
        primary_embeds = self.primary_embedding(primary)
        primary_embeds = primary_embeds.permute(0, 2, 1)
        return primary_embeds

class EdgeNetwork(nn.Module):
    def __init__(self, in_node_features, in_edge_features, out_features):
        super(EdgeNetwork, self).__init__()
        self.in_node_features = in_node_features
        self.in_edge_features = in_edge_features
        self.out_features = out_features

        # all the features are not used in every model
        # depending on the initialization, different set/combination of features are used in each model
        # we build an ensemble of models and expect different combinations to be used
        self.fc1_ef = nn.Linear(in_features=in_edge_features, out_features=8)

        self.fc1 = nn.Linear(in_features=in_node_features + 8, out_features=12)
        self.elu1 = nn.ELU()

        self.fc2 = nn.Linear(in_features=12, out_features=8)
        self.elu2 = nn.ELU()

        self.fc3 = nn.Linear(in_features=8, out_features=out_features, bias=False)

    def forward(self, x):
        ef = x[:, -self.in_edge_features:]
        ef = self.fc1_ef(ef)

        x = torch.cat((x[:, :-self.in_edge_features], ef), dim=1)

        x = self.fc1(x)
        x = self.elu1(x)
        x = self.fc2(x)
        x = self.elu2(x)
        x = self.fc3(x)
        return x

class SCONESNetwork(BaseNetwork):
    def __init__(self, *args, **kwargs):
        BaseNetwork.__init__(self, *args, **kwargs)

        self.embed_size = 6
        self.num_aa_props = 1

        self.num_node_features = self.embed_size + self.num_aa_props
        self.num_edge_features = 16

        self.embeddingNet = EmbeddingNetwork(self.embed_size)       
        self.edgeNet = EdgeNetwork(self.num_node_features * 2, self.num_edge_features, 1)
 
    def forward_single(self, primary_rt, primary_props_rt, primary_mt, primary_props_mt, neighbors, position, edge_features, edge_features_rev):
        num_neighbors = len(neighbors)
        
        # [embed_size, seqlen]
        rt_embedding = self.embeddingNet(primary_rt).reshape(self.embed_size, -1)
        mt_embedding = self.embeddingNet(primary_mt).reshape(self.embed_size, -1)

        # [num_node_features, seqlen]
        rt_embedding = torch.cat((rt_embedding, primary_props_rt), axis=0)
        mt_embedding = torch.cat((mt_embedding, primary_props_mt), axis=0)

        # [1, num_node_features]
        central_rt = rt_embedding[:, position]    
        central_mt = mt_embedding[:, position]

        # [num_neighbors, num_node_features]
        central_rt_expanded = central_rt.expand(num_neighbors, -1)
        central_mt_expanded = central_mt.expand(num_neighbors, -1)

        # [num_neighbors, num_node_features]
        rt_neighbors_embeds = torch.transpose(rt_embedding[:, neighbors], 0, 1)
        mt_neighbors_embeds = torch.transpose(mt_embedding[:, neighbors], 0, 1)

        # [num_neighbors, 2 * num_node_features + num_edge_features]
        rt_edges = torch.cat((central_rt_expanded, rt_neighbors_embeds, edge_features), axis=1)
        rt_edges_reverse = torch.cat((rt_neighbors_embeds, central_rt_expanded, edge_features_rev), axis=1)
        mt_edges = torch.cat((central_mt_expanded, mt_neighbors_embeds, edge_features), axis=1)
        mt_edges_reverse = torch.cat((mt_neighbors_embeds, central_mt_expanded, edge_features_rev), axis=1)

        # merge all vectors so that we can compute all contributions in one big batch
        # [2 * 2 * num_neighbors, a + b + c]        
        edges = torch.cat((rt_edges, rt_edges_reverse, mt_edges, mt_edges_reverse), axis=0)

        # compute edge reduction
        edges_out = self.edgeNet(edges)
        edges_out = edges_out.reshape(2, 2 * num_neighbors, -1)

        rt_output = torch.sum(edges_out[0], axis=0) * 0.25
        mt_output = torch.sum(edges_out[1], axis=0) * 0.25

        edges_out = edges_out.reshape(2, 2, num_neighbors, -1)
        equality_loss = torch.sum((edges_out[:, 0, :, :] - edges_out[:, 1, :, :])**2)

        edges_out = torch.mean(edges_out, dim=1).reshape(2, num_neighbors, -1)

        net_additive_loss = equality_loss

        ddG = (mt_output - rt_output).reshape(1, 1)
        return edges_out, net_additive_loss, ddG

    def forward(self, primary_rt, primary_props_rt, primary_mt, primary_props_mt, position, neighbors, edge_features, edge_features_rev):
        batch_size = len(primary_rt)

        ddG = torch.zeros((batch_size, 1), device=self.device)
        additive_loss = torch.zeros((batch_size, 1), device=self.device)
        edges_out = [None] * batch_size
        for b in range(batch_size):
            edges_out[b], additive_loss[b], ddG[b, 0] = self.forward_single(primary_rt[b], primary_props_rt[b], primary_mt[b], primary_props_mt[b], neighbors[b], position[b], edge_features[b], edge_features_rev[b])
        return edges_out, additive_loss, ddG
