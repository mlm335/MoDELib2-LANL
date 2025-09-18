import numpy as np
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
from matplotlib import pyplot as plt
from matplotlib import cm
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
import layer

class GCN_Inv_Phys(nn.Module):
    def __init__(self, in_feats_GP, in_feats_Node, latent_feat, out_feats_GP, out_feats_Node):
        super(GCN_Inv_Phys, self).__init__()
        self.conv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node, latent_feat, latent_feat)
        self.conv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.linear = nn.Linear(latent_feat, out_feats_Node)

        self.Bconv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node+1, latent_feat, latent_feat)
        self.Bconv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.Blinear = nn.Linear(latent_feat, 1)

        self.Econv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node, latent_feat, latent_feat)
        self.Econv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.Elinear = nn.Linear(latent_feat, 1)

        
        self.sigconv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node+1, latent_feat, latent_feat)
        self.sigconv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.siglinear = nn.Linear(latent_feat, 1)
        
    def forward(self, g, feature_GP, feature_Node, feature_tan, feature_T):
        h1,h2 = self.conv1(g, feature_GP, feature_Node)
        h1,h2 = self.conv2(g, h1,h2)
        coeff = self.linear(h2)
        out = coeff[:,0:1]*feature_tan[:,0]+ coeff[:,1:]*feature_tan[:,1]
        
        h1,h2 = self.Econv1(g, feature_GP, feature_Node)
        h1,h2 = self.Econv2(g, h1,h2)
        E = th.abs(self.Elinear(h2))
        E = th.exp(-E/feature_T)

        featT = th.cat([feature_Node, feature_T], axis=1)
        h1,h2 = self.Bconv1(g, feature_GP, featT)
        h1,h2 = self.Bconv2(g, h1,h2)
        B = th.abs(self.Blinear(h2))

          
        sigh1,sigh2 = self.sigconv1(g, feature_GP, featT)
        sigh1,sigh2 = self.sigconv2(g, sigh1,sigh2)
        sigout = self.siglinear(sigh2)
        return out/B*E, th.abs(sigout)

class GCN_Inv_Phys2(nn.Module):
    def __init__(self, in_feats_GP, in_feats_Node, latent_feat, out_feats_GP, out_feats_Node):
        super(GCN_Inv_Phys2, self).__init__()
        self.conv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node, latent_feat, latent_feat)
        self.conv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.linear = nn.Linear(latent_feat, out_feats_Node)

        self.Bconv1 = layer.GCNLayer_DDD(1, in_feats_Node+1, latent_feat, latent_feat)
        self.Bconv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.Blinear = nn.Linear(latent_feat, 1)

        self.Econv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node, latent_feat, latent_feat)
        self.Econv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.Elinear = nn.Linear(latent_feat, 1)

        
        self.sigconv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node+1, latent_feat, latent_feat)
        self.sigconv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.siglinear = nn.Linear(latent_feat, 1)
        
    def forward(self, g, feature_GP, feature_Node,feature_bdott, feature_tan, feature_T):
        h1,h2 = self.conv1(g, feature_GP, feature_Node)
        h1,h2 = self.conv2(g, h1,h2)
        coeff = self.linear(h2)
        out = coeff[:,0:1]*feature_tan[:,0]+ coeff[:,1:]*feature_tan[:,1]
        
        h1,h2 = self.Econv1(g, feature_GP, feature_Node)
        h1,h2 = self.Econv2(g, h1,h2)
        E = th.abs(self.Elinear(h2))
        E = th.exp(-E/feature_T)

        featT = th.cat([feature_Node, feature_T], axis=1)
        h1,h2 = self.Bconv1(g, feature_bdott, featT)
        h1,h2 = self.Bconv2(g, h1,h2)
        B = th.abs(self.Blinear(h2))

          
        sigh1,sigh2 = self.sigconv1(g, feature_GP, featT)
        sigh1,sigh2 = self.sigconv2(g, sigh1,sigh2)
        sigout = self.siglinear(sigh2)
        return out/B*E, th.abs(sigout)
    def BEV(self, g, feature_GP, feature_Node,feature_bdott, feature_tan, feature_T):
        h1,h2 = self.conv1(g, feature_GP, feature_Node)
        h1,h2 = self.conv2(g, h1,h2)
        coeff = self.linear(h2)
        out = coeff[:,0:1]*feature_tan[:,0]+ coeff[:,1:]*feature_tan[:,1]
        
        h1,h2 = self.Econv1(g, feature_GP, feature_Node)
        h1,h2 = self.Econv2(g, h1,h2)
        E = th.abs(self.Elinear(h2))
#        E = th.exp(-E/feature_T)

        featT = th.cat([feature_Node, feature_T], axis=1)
        h1,h2 = self.Bconv1(g, feature_bdott, featT)
        h1,h2 = self.Bconv2(g, h1,h2)
        B = th.abs(self.Blinear(h2))

          
        sigh1,sigh2 = self.sigconv1(g, feature_GP, featT)
        sigh1,sigh2 = self.sigconv2(g, sigh1,sigh2)
        sigout = self.siglinear(sigh2)
        return out, B, E
class GCN_Inv_Phys2_MDLoopDD(nn.Module):
    def __init__(self, in_feats_GP, in_feats_Node, latent_feat, out_feats_GP, out_feats_Node):
        super(GCN_Inv_Phys2_MDLoopDD, self).__init__()
        self.conv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node, latent_feat, latent_feat)
        self.conv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.linear = nn.Linear(latent_feat, out_feats_Node)

        self.Bconv1 = layer.GCNLayer_DDD(1, in_feats_Node+1, latent_feat, latent_feat)
        self.Bconv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.Blinear = nn.Linear(latent_feat, 1)

        self.Econv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node, latent_feat, latent_feat)
        self.Econv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.Elinear = nn.Linear(latent_feat, 1)

        
        self.sigconv1 = layer.GCNLayer_DDD(in_feats_GP, in_feats_Node+1, latent_feat, latent_feat)
        self.sigconv2 = layer.GCNLayer_DDD(latent_feat, latent_feat, latent_feat, latent_feat)
        self.siglinear = nn.Linear(latent_feat, 1)
        
    def forward(self, g, feature_GP, feature_Node,feature_bdott, feature_tan, feature_T):
        h1,h2 = self.conv1(g, feature_GP, feature_Node)
        h1,h2 = self.conv2(g, h1,h2)
        coeff = self.linear(h2)
        out = coeff[:,0:1]*feature_tan[:,0]+ coeff[:,1:2]*feature_tan[:,1]+ coeff[:,2:3]*feature_tan[:,2]
        
        h1,h2 = self.Econv1(g, feature_GP, feature_Node)
        h1,h2 = self.Econv2(g, h1,h2)
        E = th.abs(self.Elinear(h2))
        E = th.exp(-E/feature_T)

        featT = th.cat([feature_Node, feature_T], axis=1)
        h1,h2 = self.Bconv1(g, feature_bdott, featT)
        h1,h2 = self.Bconv2(g, h1,h2)
        B = th.abs(self.Blinear(h2))

          
        sigh1,sigh2 = self.sigconv1(g, feature_GP, featT)
        sigh1,sigh2 = self.sigconv2(g, sigh1,sigh2)
        sigout = self.siglinear(sigh2)
        return out/B*E, th.abs(sigout)
    def BEV(self, g, feature_GP, feature_Node,feature_bdott, feature_tan, feature_T):
        h1,h2 = self.conv1(g, feature_GP, feature_Node)
        h1,h2 = self.conv2(g, h1,h2)
        coeff = self.linear(h2)
        out = coeff[:,0:1]*feature_tan[:,0]+ coeff[:,1:2]*feature_tan[:,1]+ coeff[:,2:3]*feature_tan[:,2]
        
        h1,h2 = self.Econv1(g, feature_GP, feature_Node)
        h1,h2 = self.Econv2(g, h1,h2)
        E = th.abs(self.Elinear(h2))
#        E = th.exp(-E/feature_T)

        featT = th.cat([feature_Node, feature_T], axis=1)
        h1,h2 = self.Bconv1(g, feature_bdott, featT)
        h1,h2 = self.Bconv2(g, h1,h2)
        B = th.abs(self.Blinear(h2))

          
        sigh1,sigh2 = self.sigconv1(g, feature_GP, featT)
        sigh1,sigh2 = self.sigconv2(g, sigh1,sigh2)
        sigout = self.siglinear(sigh2)
        return out, B, E

