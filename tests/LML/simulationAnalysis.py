# /opt/local/bin/python3.12 test.py
import sys, string, os, re, glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import matplotlib as mpl
from collections import defaultdict
import matplotlib.cm as cm
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
sys.path.append("../../python/lib")
from readFile import * 
from readF import *
from modlibUtils import *
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib
from utilities import *


# -----------------------------
# User Settings
# -----------------------------
# simulationDir = os.path.join("../..", "tutorials", "LML_MD_Copy_LargeBox_CS1_mu49_Lmin5_Lmax15_QP0.05_RM50")  # path to your simulation directory, where inputFiles/ and evl/ are located
simulationDir = "/Users/matthewmaron/Documents/MoDELib2-LANL/tutorials/LML_MD_Copy_LargerBox_CS1_mu49_Lmin5_Lmax15_QP0.05_RM20"
simulation_type = "compression" # "tension" or "compression"

# --- Segment Color Scheme ---
colorScheme = "characterAngle" # options: "characterAngle"/ "glissileSessile"
# Note: sessile segments always colored as magenta 
# Note: glissileSessile will color glissile as blue, sessile as magenta
colorBar = False
cmap = plt.cm.viridis
norm = mpl.colors.Normalize(vmin=0.0, vmax=180.0)
# -----------------------------
# Direcctories 
# -----------------------------
evl_dir = os.path.join(simulationDir, "evl")
out_root = os.path.join(os.getcwd(), os.path.basename(simulationDir))  # "./LML_MD_Copy"
os.makedirs(out_root, exist_ok=True)
fig_dir  = os.path.join(out_root, "figures")
os.makedirs(fig_dir, exist_ok=True)
out_txt = os.path.join(out_root, "simulation_summary.txt")

# -----------------------------
# Setup MoDELib objects once 
# -----------------------------
ddBase = pyMoDELib.DislocationDynamicsBase(simulationDir)
ddio   = pyMoDELib.DDconfigIO(evl_dir)
defectiveCrystal = pyMoDELib.DefectiveCrystal(ddBase)
evl_steps = list_evl_steps(evl_dir)
xMin = np.array(ddBase.mesh.xMin(), dtype=float)
xMax = np.array(ddBase.mesh.xMax(), dtype=float)

# -----------------------------
# Read matfile once
# -----------------------------
matfile = os.path.join(simulationDir, 'inputFiles', 'Fe_320.txt')
mu0_SI = get_scalar(matfile, 'mu0_SI')
rho_SI = get_scalar(matfile, 'rho_SI')
b_SI = get_scalar(matfile, 'b_SI')
v_dd2SI = np.sqrt(mu0_SI / rho_SI)
t_dd2SI = b_SI / v_dd2SI

# -----------------------------
# Read polyfile once
# -----------------------------
polyfile = os.path.join(simulationDir, 'inputFiles', 'polycrystal.txt')
deformationGradient=get_matrix(polyfile,'F')
V_b3 = abs(np.linalg.det(deformationGradient))   # volume in b^3
densityFactor = 1.0 / (V_b3 * b_SI**2)

# -----------------------------
# Read F once
# -----------------------------
F, Flabels = readFfile(simulationDir)
runID = getFarray(F, Flabels, 'runID')
t = getFarray(F, Flabels, 'time [b/cs]')
e_33 = getFarray(F, Flabels, 'e_33')
s_33 = getFarray(F, Flabels, 's_33')
strain = np.array(e_33, dtype=float)          # dimensionless
stress = np.array(s_33, dtype=float)          # whatever units your F-file uses

# -----------------------------
# Map EVL -> second index in F-file arrays
# -----------------------------
runID_arr = np.asarray(runID, dtype=int)
# collect all indices per runID
rid_to_indices = defaultdict(list)
for i, rid in enumerate(runID_arr):
    rid_to_indices[int(rid)].append(i)
def second_index_for_runid(evl):
    idxs = rid_to_indices.get(int(evl), [])
    if len(idxs) == 0:
        return None
    if len(idxs) == 1:
        return idxs[0]          # only one available
    return idxs[1]              # ALWAYS the second occurrence


# -----------------------------
# Load configuration
# -----------------------------
with open(out_txt, "w") as f_out:
    for evl_file in evl_steps:

        # --- match to F-file row (second occurrence if duplicated) ---
        iF = second_index_for_runid(evl_file)
        if iF is None:
            print(f"[warn] EVL {evl_file}: not found in F-file runID, skipping row.")
            continue
        time_s     = float(t[iF]) * t_dd2SI
        strain_pct = float(strain[iF]) * 100.0
        stress_GPa = float(stress[iF]) * mu0_SI * 1e-9
        if simulation_type == "compression":
            strain_pct = -strain_pct
            stress_GPa = -stress_GPa

        # --- load EVL config ---
        ddio.readTxt(evl_file)
        defectiveCrystal.initializeConfiguration(ddio)
        DN = defectiveCrystal.dislocationNetwork()

        # --- segments + lengths ---
        segments = dislocationSegments_from_Network(DN)
        Lgl, Lse, Lbnd, Lgb, Lscrew, Ledge, Lmixed = network_length_from_networklinks(DN,screw_dev=10.0,edge_lo=80.0,edge_hi=100.0)
        Ltot = Lgl + Lse + Lbnd + Lgb

        # --- write summary info ---
        # --- densities in SI (m^-2), using your densityFactor = 1/(V_b3 * b_SI^2) ---
        rho_tot    = Ltot   * densityFactor
        rho_gl     = Lgl    * densityFactor
        rho_se     = Lse    * densityFactor
        rho_bnd    = Lbnd   * densityFactor
        rho_gb     = Lgb    * densityFactor
        rho_screw  = Lscrew * densityFactor
        rho_edge   = Ledge  * densityFactor
        rho_mixed  = Lmixed * densityFactor
        # --- write one line (no header) ---
        f_out.write(
            f"{evl_file:d} "
            f"{time_s:.12e} "
            f"{strain_pct:.8f} "
            f"{stress_GPa:.12e} "
            f"{rho_tot:.12e} "
            f"{rho_gl:.12e} "
            f"{rho_se:.12e} "
            f"{rho_bnd:.12e} "
            f"{rho_gb:.12e} "
            f"{rho_screw:.12e} "
            f"{rho_edge:.12e} "
            f"{rho_mixed:.12e}\n"
        )

        # --- plot Segments and Save Figure ---
        fig = plt.figure(figsize=(6, 12))
        ax = fig.add_subplot(111, projection="3d")
        for (lk, a, b, type, theta) in segments:
            if colorScheme == "characterAngle":
                if type == "sessile" or theta is None:
                    color = "m"       # sessile / undefined
                    alpha = 0.9
                else:
                    theta = float(np.clip(theta, 0.0, 180.0))
                    color = cmap(norm(theta))
                    alpha = 1.0
                ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]],lw=1.5, alpha=alpha, color=color)

            elif colorScheme == "glissileSessile":
                if type=="glissile":
                    ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], lw=1.5, alpha=1.0, color='b')
                else:
                    ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], lw=1.5, alpha=1.0, color='m')

        draw_box(ax, xMin, xMax, lw=1.5, alpha=0.5)
        ax.set_box_aspect((xMax - xMin))
        ax.set_axis_off()

        if colorScheme == "characterAngle" and colorBar:
            sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            # cbar = plt.colorbar(sm, ax=ax, pad=0.08, shrink=0.75)
            cbar = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02, shrink=0.92, aspect=35)
            cbar.set_label("Character angle (deg)")
            cbar.set_ticks([0, 45, 90, 135, 180])
            sessile_handle = Line2D([0],[0], color='m', lw=2, label='sessile')
            ax.legend(handles=[sessile_handle], loc='upper left', framealpha=0.85)


        # ax.set_title(rf"EVL {evl_file}")
        fig_path = os.path.join(fig_dir, f"evl_{evl_file:08d}.png")
        plt.savefig(fig_path, dpi=200, bbox_inches="tight")
        plt.close(fig)


print(f"Wrote summary: {out_txt}")
print(f"Saved figures to: {fig_dir}")

