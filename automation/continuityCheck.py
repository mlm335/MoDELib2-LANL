import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
sys.path.append("../python")
from modlibUtils import *
from stereoUtils import *
sys.path.append("../build/tools/pyMoDELib")
import pyMoDELib

# ----------------------------------------------------------------- #
def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def compute_rotation_matrix(x1_norm,x3_norm):
    x3_dir = x3_norm
    x1_dir = x1_norm - np.dot(x1_norm, x3_dir) * x3_dir # orthogonalize x1 to x3
    x2_dir = np.cross(x3_dir, x1_dir)
    return np.row_stack((x1_dir, x2_dir, x3_dir))

def stress_to_voigt(sigma):
    return [ sigma[0,0], sigma[1,1], sigma[2,2], sigma[0,1], sigma[1,2], sigma[0,2] ]
# ----------------------------------------------------------------- #

# -- Material File Info -- #
materialFile="../Library/Materials/Fe_320.txt"
b_SI=getValueInFile(materialFile,'b_SI')
mu_SI=getValueInFile(materialFile,'mu0_SI')
rho_SI=getValueInFile(materialFile,'rho_SI')
cs = np.sqrt(mu_SI/rho_SI)  

# -- Inputs to simulations -- #
X1_Direction = [ [1,1,-1],[2,3,-2],[1,2,-1],[1,3,-1],[1,4,-1],[0,1,0],[-1,2,1],[-1,1,1],[-2,1,2],[-1,0,1],[-2,1,2],[-3,-2,3] ]
# Resolved_Shear = [0.0122, 0.0091463, 0.0061, 0.0031, 0.00122]
Temperature = [100, 300, 500, 700, 1000]
Resolved_Shear = np.linspace(0.01, 1, 1000)/mu_SI/1e-9

# -- Local Reference Frame Unit Vectors -- #
b_unit = normalize([1,1,-1])
n_unit = normalize([1,0,1])

# -- Begin Check up on Mobility -- #
results = {}
angles_list = []
results = {temp: [] for temp in Temperature}
for x1 in X1_Direction:
    x1_unit = normalize(x1)
    cos_theta = np.clip(np.dot(b_unit, x1_unit), -1.0, 1.0)
    miss_angle = float(np.degrees(np.arccos(cos_theta)))
    C2G = compute_rotation_matrix(x1_unit, n_unit)  
    b_g = C2G @ b_unit
    x1_g= C2G @ x1_unit
    n_g = C2G @ n_unit
    
    results = {temp: [] for temp in Temperature}
    for T in Temperature:
        mat=pyMoDELib.PolycrystallineMaterialBase(materialFile,T) 
        DDPy="../python/MobilityLaw/DDMobilityPy" 
        mobPy=pyMoDELib.DislocationMobilityPy(mat,DDPy) 

        for s0 in Resolved_Shear:
            S_crystal = s0* (np.outer(b_unit, n_unit)+np.outer(n_unit, b_unit))
            S_global = C2G @ S_crystal @ C2G.T
            v = mobPy.velocity(S_global, b_g, x1_g, n_g, T) # Code Units 
            results[T].append(v*cs) # convert to m/s 

    # -- Plot
    cmap = plt.get_cmap("viridis")
    color_map = {temp: cmap(i / len(Temperature)) for i, temp in enumerate(Temperature)}
    plt.figure(figsize=(8, 8))
    for temp, velocities in results.items():
        plt.plot(Resolved_Shear*mu_SI*1e-9, velocities, color=color_map[temp], label=f'{temp}K')
    plt.xlabel('Resolved Shear Stress [GPa]',fontsize=16)
    plt.ylabel('Predicted Velocity [A/ps]',fontsize=16)
    plt.title(rf'$\theta={miss_angle:.1f}^\circ$')
    plt.tick_params(axis = 'both', which = 'major', direction='in', bottom=True, left=True, top=True, right=True, labelsize = 14)
    plt.tick_params(axis = 'both', which = 'minor', direction='in', bottom=True, left=True, top=True, right=True, labelsize = 14)
    plt.legend()
    plt.grid()
    plt.savefig(f"ContinuityPlots/Velocity_X1_{x1[0]}{x1[1]}{x1[2]}.pdf", dpi=1000)
# plt.show()



