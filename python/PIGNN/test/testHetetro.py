import torch as th
import dgl
import dgl.function as fn
import torch.nn as nn
import torch.nn.functional as F
import numpy as np



class GCNLayer_DDD(nn.Module):
    def __init__(self, in_feats_GP, in_feats_Node, out_feats_GP, out_feats_Node):
        super(GCNLayer_DDD, self).__init__()
        self.linear_Node = nn.Linear(in_feats_GP+in_feats_Node, out_feats_Node)
        self.linear_GP = nn.Linear(out_feats_Node, out_feats_GP)
    def forward(self, g, feature_GP, feature_Node):
        # Creating a local scope so that all the stored ndata and edata                                                                             
        # (such as the `'h'` ndata below) are automatically popped out                                                                              
        # when the scope exits.                                                                                                                     
        with g.local_scope():
            g.nodes['GP'].data['h'] = feature_GP
            g.nodes['Node'].data['h'] = feature_Node
            g.update_all(fn.copy_u('h','m'), fn.mean('m','hGP'), etype ='edge')
            h = th.cat((g.nodes['Node'].data['h'],g.nodes['Node'].data['hGP']),axis=1)

            g.nodes['Node'].data['h'] = F.relu(self.linear_Node(h))
            g.update_all(fn.copy_u('h','m'), fn.mean('m','hNode'), etype = 'iedge')
            h = g.nodes['GP'].data['hNode']
            g.nodes['GP'].data['h'] = F.relu(self.linear_GP(h))
            return g.nodes['GP'].data['h'], g.nodes['Node'].data['h']

class test_Hetero():
    def __init__(self):
        super().__init__()
##### construct a simple heterogeneous graph #####

        dst = np.array([0,1,1,2])
        src = np.array([0,0,1,1])
        
        isrc = np.array([0,1,1,2])
        idst = np.array([0,0,1,1])
        features_GP = np.random.randn(2,7)
        features_Node = np.random.randn(3,3)
        
        self.graph = dgl.heterograph({('GP','edge','Node'):(src,dst), ('Node','iedge','GP'):(isrc,idst)})
        self.graph.nodes['GP'].data['GPfeat'] = th.from_numpy(features_GP)
        self.graph.nodes['Node'].data['Nodefeat'] = th.from_numpy(features_Node)
#### Single-layer HeteroSage ####
        self.model  = GCNLayer_DDD(7,3,2,2).double()

        
    def forward(self):
        with th.no_grad():
            features_GP = self.graph.nodes['GP'].data['GPfeat']
        
            features_node = self.graph.nodes['Node'].data['Nodefeat']

            h1, h2 = self.model(self.graph, features_GP, features_node)

        return h1.numpy()

