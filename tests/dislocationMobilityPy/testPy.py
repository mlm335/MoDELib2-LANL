# /opt/local/bin/python3.13 test.py
import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
from stereoUtils import *
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib


def construct_stress_tensor_general(resolved_shear, b_unit, n_unit):
    S = np.outer(b_unit, n_unit) + np.outer(n_unit, b_unit)  # symmetric shear tensor
    S /= np.linalg.norm(S)  # normalize to unit Frobenius norm
    return resolved_shear * S

def compute_rotation_matrix(x1_norm,x3_norm):
    x3_dir = x3_norm
    x1_dir = x1_norm - np.dot(x1_norm, x3_dir) * x3_dir # orthogonalize x1 to x3
    x2_dir = np.cross(x3_dir, x1_dir)
    return np.row_stack((x1_dir, x2_dir, x3_dir))


materialFile="../../Library/Materials/Fe_320.txt" # material file
b_SI=getValueInFile(materialFile,'b_SI')
mu_SI=getValueInFile(materialFile,'mu0_SI')
rho_SI=getValueInFile(materialFile,'rho_SI')
cs = np.sqrt(mu_SI/rho_SI)  # speed of sound in material
T_array = [100, 300, 500, 700, 1000]  # K
s0_array = np.linspace(0.01, 1, 1000)/mu_SI/1e-9

x1=np.array([1,1,1]) # x1 crystal direction
x1=x1/np.linalg.norm(x1); 
b=np.array([1,1,1]) # burgers vector
b=b/np.linalg.norm(b)
n=np.array([-1,0,1]) # plane normal
n=n/np.linalg.norm(n)

results = {temp: [] for temp in T_array}
for T in T_array:
    mat=pyMoDELib.PolycrystallineMaterialBase(materialFile,T) # material object
    DDPy="../../python/MobilityLaw/DDMobilityPy" 
    mobPy=pyMoDELib.DislocationMobilityPy(mat,DDPy) # mobility object
    for s0 in s0_array:
        S_crystal = construct_stress_tensor_general(s0,b,n)
        C2G = compute_rotation_matrix(x1, n)  # rotation matrix from crystal to global coordinates
        S_global = C2G @ S_crystal @ C2G.T
        # v = mobPy.velocity(S_crystal, b, C2G[1,:], n, T)
        v = mobPy.velocity(S_global, C2G @ b, C2G @ x1, C2G @ n, T)
        results[T].append(v*cs/100) # convert to A/ps 


# -- Plot
cmap = plt.get_cmap("viridis")
color_map = {temp: cmap(i / len(T_array)) for i, temp in enumerate(T_array)}
plt.figure(figsize=(8, 8))
for temp, velocities in results.items():
    plt.plot(s0_array*mu_SI*1e-9, velocities, color=color_map[temp], label=f'{temp}K')
plt.xlabel('Resolved Shear Stress [GPa]',fontsize=16)
plt.ylabel('Predicted Velocity [A/ps]',fontsize=16)
plt.tick_params(axis = 'both', which = 'major', direction='in', bottom=True, left=True, top=True, right=True, labelsize = 14)
plt.tick_params(axis = 'both', which = 'minor', direction='in', bottom=True, left=True, top=True, right=True, labelsize = 14)
plt.legend()
plt.grid()
plt.show()


# #### Dr. Po's Test ####
# m=np.array([0,0,1]) # stress axis
# m=m/np.linalg.norm(m); # normalized stress axis
# s0=0.01 # stress amplitude
# S=np.outer(m,m)*s0 # stress tensor
# b=np.array([1,1,1]) # burgers vector
# b=b/np.linalg.norm(b)
# n=np.array([-1,0,1]) # plane normal
# n=n/np.linalg.norm(n)
# for thetaDeg in range(90+1):
# # thetaDegrees=[0,90]
# # for thetaDeg in thetaDegrees:
#     theta=thetaDeg*np.pi/180 # angle between line tangent and burgers vector
#     xi=angleAxis(theta,n)@b # rotate b about n to define the line tangent
#     tau=n@S@b.transpose() # shear stress
#     # print(S, b, xi, n, T)
#     v=mobPy.velocity(S,b,xi,n,T)
#     print(v)

