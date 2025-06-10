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

materialFile="../../Library/Materials/W.txt"
T=300.0
mat=pyMoDELib.PolycrystallineMaterialBase(materialFile,T) # material object
mob111=pyMoDELib.DislocationMobilityBCC(mat) # mobility object

m=np.array([0,0,1]) # stress axis
m=m/np.linalg.norm(m); # normalized stress axis
s0=0.01 # stress amplitude
S=np.outer(m,m)*s0 # stress tensor
b=np.array([1,1,1]) # burgers vector
b=b/np.linalg.norm(b)
n=np.array([-1,0,1]) # plane normal
n=n/np.linalg.norm(n)
for thetaDeg in range(90+1):
    theta=thetaDeg*np.pi/180 # angle between line tangent and burgers vector
    xi=angleAxis(theta,n)@b # rotate b about n to define the line tangent
    tau=n@S@b.transpose() # shear stress
    v=mob111.velocity(S,b,xi,n,T)
    print(v)

# STARTING TO MAKE stereographic triangle plot

# Corners of the stereographic triangle
d1=np.array([0,1,1]);
d1=d1/np.linalg.norm(d1);
d2=np.array([-1,1,1]);
d2=d2/np.linalg.norm(d2);
d3=np.array([0,0,1]);
d3=d3/np.linalg.norm(d3);

# Stereographic grid
numPts=20
theta=np.linspace(0,numPts,numPts,endpoint=True)/numPts*np.acos(np.dot(d2,d3));
phi  =np.linspace(0,numPts,numPts,endpoint=True)/numPts*np.pi/4;
[Theta,Phi]=np.meshgrid(theta,phi);
Z=np.cos(Theta);
X=-np.sin(Theta)*np.sin(Phi);
Y=np.sin(Theta)*np.cos(Phi);
C = cm.jet(np.random.rand(Z.shape[0],Z.shape[1]))
C=cm.jet(Z)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#surf = ax.plot_surface(X, Y, Z,Z,linewidth=0, antialiased=False)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors = C,
        linewidth=0, antialiased=False)
plt.show()
