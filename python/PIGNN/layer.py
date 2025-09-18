import numpy as np
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
import glob
from time import time

import dgl
from dgl.nn import GraphConv
import dgl.function as fn
import networkx as nx

import torch as th
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd.functional import vjp, vhp, jacobian, hessian
# weird MKL issue started 8/16 after restart...
import os

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
class GCNLayer_GP(nn.Module):
    def __init__(self, in_feats_GP,  out_feats_GP, out_feats_Node):
        super(GCNLayer_GP, self).__init__()
        self.linear_Node = nn.Linear(in_feats_GP, out_feats_Node)
        self.linear_GP = nn.Linear(out_feats_Node, out_feats_GP)
    def forward(self, g, feature_GP):
        # Creating a local scope so that all the stored ndata and edata
        # (such as the `'h'` ndata below) are automatically popped out
        # when the scope exits.
        with g.local_scope():
            g.nodes['GP'].data['h'] = feature_GP

            g.update_all(fn.copy_u('h','m'), fn.mean('m','hGP'), etype ='edge')
            h = g.nodes['Node'].data['hGP']
            g.nodes['Node'].data['h'] = F.relu(self.linear_Node(h))
            g.update_all(fn.copy_u('h','m'), fn.mean('m','hNode'), etype = 'iedge')
            h = g.nodes['GP'].data['hNode']
            g.nodes['GP'].data['h'] = F.relu(self.linear_GP(h))
            return g.nodes['GP'].data['h'], g.nodes['Node'].data['h']
        
