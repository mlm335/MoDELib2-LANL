# /opt/local/bin/python3.13 test.py
import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib

materialFile="../../Library/Materials/Zr.txt"
absoluteTemperature=0.0
mat=pyMoDELib.PolycrystallineMaterialBase(materialFile,absoluteTemperature)
energyDensityScaling=mat.mu_SI*mat.b_SI
C2G=np.array([[1,0,0],[0,1,0],[0,0,1]])
singleCrystal=pyMoDELib.SingleCrystalBase(mat,C2G)

# initialize
n0=2 # image repetitions along first direction
n1=2 # image repetitions along second direction
pts=50 # GSFE resolution in each direction
X=np.zeros((pts+1,pts+1))
Y=np.zeros((pts+1,pts+1))
GSFE=dict() # stores GSFE of matrix and other phases
GSFE["matrix"]=np.zeros((pts+1,pts+1))
for phaseKey, phaseValue in singleCrystal.secondPhases().items():
   GSFE[phaseValue.name]=np.zeros((pts+1,pts+1))

# Compute and plot
for planeKey, planeValue in singleCrystal.planeNormals().items(): # loop over lattice planes
    v0=planeValue.primitiveVectors[0].cartesian() # first primitive vector of the matrix lattice plane
    v1=planeValue.primitiveVectors[1].cartesian() # second primitive vector of the matrix lattice plane
    for i in range(pts+1):
        for j in range(pts+1):
            s=i*n0/pts*v0+j*n1/pts*v1 # slip vector
            sL=planeValue.localSlipVector(s) # slip vector in the plane local coordinates
            X[j,i]=sL[0]
            Y[j,i]=sL[1]
            GSFE["matrix"][j,i]=planeValue.misfitEnergy(s)*energyDensityScaling
            for phaseKey, phaseValue in singleCrystal.secondPhases().items():
                GSFE[phaseValue.name][j,i]=phaseValue.misfitEnergy(s,planeKey)*energyDensityScaling

    # Plot GSFE surface
    fig, ax = plt.subplots()
    plt.gca().set_aspect("equal")
    pmsh=ax.pcolormesh (X, Y, GSFE["matrix"],cmap=cm.coolwarm,shading='gouraud')
#    fig.colorbar(pmsh)
    fig.savefig(mat.materialName+"_plane_"+str(planeKey)+"_matrix", bbox_inches='tight')
    for phaseKey, phaseValue in singleCrystal.secondPhases().items():
        pmsh=ax.pcolormesh (X, Y, GSFE[phaseValue.name],cmap=cm.coolwarm,shading='gouraud')
#        fig.colorbar(pmsh)
        fig.savefig(mat.materialName+"_plane_"+str(planeKey)+"_"+phaseValue.name, bbox_inches='tight')
