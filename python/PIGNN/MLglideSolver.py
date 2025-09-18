import torch as th
import dgl
import dgl.function as fn
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from gmodel import GCN_Inv_Phys2_MDLoopDD
from utility import GraphConstruct, VTKparser


class PIGNN():
    def __init__(self):
        super().__init__()
##### construct a simple heterogeneous graph #####

        dst = np.array([0,1,1,2])
        src = np.array([0,0,1,1])
        
        isrc = np.array([0,1,1,2])
        idst = np.array([0,0,1,1])
        features_GP = np.random.randn(2,7)
        features_Node = np.random.randn(3,3)
        features_bdott = np.random.randn(2,1)
        features_tan = np.random.randn(3,3,3)
        features_T = np.random.randn(3,1)
        
        self.graph = dgl.heterograph({('GP','edge','Node'):(src,dst), ('Node','iedge','GP'):(isrc,idst)})
        self.graph.nodes['GP'].data['GPfeat'] = th.from_numpy(features_GP)
        self.graph.nodes['Node'].data['Nodefeat'] = th.from_numpy(features_Node)

        self.graph.nodes['GP'].data['bdott'] = th.from_numpy(features_bdott)
        self.graph.nodes['Node'].data['Nodetan'] = th.from_numpy(features_tan)
        self.graph.nodes['Node'].data['T'] = th.from_numpy(features_T)
################################################
        self.datareader = GraphConstruct()
        
#### Multi-layer HeteroSage ####
        self.modellist  = [GCN_Inv_Phys2_MDLoopDD(7,3,20,3,3).double() for i in range(5)] # top 25% models
        self.load_parameter()
    def read_testdata(self):
        out = VTKparser('quadrature_5000.vtk')
        pos = out['points']
        stress = out['stress']
        burgers_vector_global = out['burgers_vector_global']
        CELLS = out['CELLS']
        T = 500
 #       self.graph = self.datareader.constructgraph(pos, stress, burgers_vector_global, CELLS, T)
        return self.MLmobility(pos, stress, burgers_vector_global, CELLS, T)
    def MLmobility(self, pos, stress, burgers_vector_global, CELLS, T):
        self.graph = self.datareader.constructgraph(pos, stress, burgers_vector_global, CELLS, T)
        return self.forward()

    def load_parameter(self):
        for i in range(len(self.modellist)):
            self.modellist[i].load_state_dict(th.load("../../../python/PIGNN/model_MDLoop"+str(i)+".params",map_location='cpu'))
    def forward(self):
        with th.no_grad():
            features_GP = self.graph.nodes['GP'].data['GPfeat']
        
            features_node = self.graph.nodes['Node'].data['Nodefeat']

            features_bdott = self.graph.nodes['GP'].data['bdott']

            features_tan = self.graph.nodes['Node'].data['Nodetan']

            features_T = self.graph.nodes['Node'].data['T']

            step = []
            for i in range(len(self.modellist)):
                h1, h2 = self.modellist[i](self.graph, features_GP, features_node, features_bdott, features_tan, features_T)
                step.append(h1.numpy())
        return np.stack(step).mean(axis=0)
