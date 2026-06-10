# Some useful functions
import numpy as np
import math, sys, string, os
from matplotlib import pyplot as plt
from matplotlib import rc
import pandas as pd
import matplotlib.patches as patches
from os import listdir
from os.path import isfile, join
print(os.environ['PATH'])
rc('font', family='Times New Roman')
sys.path.append('lib')
import pathlib
from readF import *
from readFile import *
from readEVL import *
from readAUX import *


def getLoopData(folderName):
    print("Getting Loop Data")
    F=np.loadtxt(folderName+'/F/F_0.txt');
    print(folderName+'/F/F_0.txt has size ' + str(np.shape(F)));
    b=[];
    n=0;
    for runID in F[:,0]:
        evl=readEVLtxt(folderName+'/evl/evl_'+str(int(runID)))
        b.append(evl.loopsBurger);
        n=n+1;
    return F[:,1], b;

def getQPData(folderName):
    print("Getting Quadrature Point Data")
    F=np.loadtxt(folderName+'/F/F_0.txt');
    print(folderName+'/F/F_0.txt has size ' + str(np.shape(F)));
    tangent=[]
    velocity=[]
    n=0;
    for runID in F[:,0]:
        aux=readAUXtxt(folderName+'/evl/ddAux_'+str(int(runID)))
        tangent.append(aux.tangent);
        velocity.append(aux.velocity);
        n=n+1;
    return F[:,1], tangent, velocity;
###########################################################################
###########################################################################
###########################################################################


parent_directory = str(pathlib.Path().resolve())+'/automationOutput'
print("Parent Directory = ", parent_directory)
dataFolders = [f for f in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, f))]

runID=0
collect_data = {}
for folder in dataFolders:
    dataFolder = os.path.join(parent_directory, folder)
    print("Data Folder = ", dataFolder)

    # -- Material properties
    mat_filepath = dataFolder + '/inputFiles/Fe_320.txt'
    print("Material File Path = ", mat_filepath)
    mu0_SI = get_scalar(mat_filepath, 'mu0_SI')
    rho_SI = get_scalar(mat_filepath, 'rho_SI')
    b_SI = get_scalar(mat_filepath, 'b_SI')
    v_dd2SI = np.sqrt(mu0_SI / rho_SI)
    t_dd2SI = b_SI / v_dd2SI
    cs = b_SI / t_dd2SI

    # -- Polycrystal
    poly_filepath = dataFolder + '/inputFiles/polycrystal.txt'
    Temperature = get_scalar(poly_filepath, 'absoluteTemperature')

    # -- ElasticDeformation
    Elastic_filepath = dataFolder + '/inputFiles/ElasticDeformation.txt'
    ExternalStress0 = get_vector(Elastic_filepath, 'ExternalStress0')
    Stress0 = np.abs(ExternalStress0[5]*mu0_SI*1e-9)

    # -- Loop and Quadrature Point Data
    t, b = getLoopData(dataFolder)
    t2, tangent, velocity = getQPData(dataFolder)

    # -- Misorientation Angle
    burgers = np.array(b[runID])
    miss_angle = []
    for k in range(len(tangent[runID])):
        QPtangent = np.array(tangent[runID][k])
        miss_product = np.dot(burgers, QPtangent)
        magnitude_a = np.linalg.norm(burgers)
        magnitude_b = np.linalg.norm(QPtangent)
        cos_theta = miss_product / (magnitude_a * magnitude_b)
        angle_radians = np.arccos(cos_theta)
        miss_angle.append(np.degrees(angle_radians))

    # -- Velocity Magnitude
    VelocityMagnitudes = np.linalg.norm(velocity[runID], axis=1)*cs # [m/s]

    # Store
    stress_str = f"{np.round(Stress0, 2)}"
    collect_data.setdefault(stress_str, []).append({
        "angles": np.array(miss_angle),
        "vels": np.array(VelocityMagnitudes),
        "temp": int(Temperature),
        "stress0": Stress0
    })


# ---------- Make one figure with 5 subplots (one per stress) ----------
stress_keys_sorted = sorted(collect_data.keys(), key=lambda s: float(s))
stress_keys_to_plot = stress_keys_sorted[:5]
n = len(stress_keys_to_plot)


# define a global color map for consistency
all_temps = sorted({int(e["temp"]) for s in stress_keys_to_plot for e in collect_data[s]})
cmap = plt.get_cmap("viridis")
color_map = {T: cmap(i / max(1, len(all_temps)-1)) for i, T in enumerate(all_temps)}

# Define custom colors for the temperatures
custom_colors = {
    100:  "tab:blue",
    300:  "tab:orange",
    500:  "tab:green",
    700:  "tab:red",
    1000: "tab:purple"
}

fig, axes = plt.subplots(1, n, figsize=(5.2*n, 4.8), sharex=True, sharey=False)
for ax, s in zip(axes, stress_keys_to_plot):

    entries = sorted(collect_data[s], key=lambda e: int(e["temp"]))
    for e in entries:
        # ax.plot(e["angles"], e["vels"], 'o', color=color_map[e["temp"]], label=f"{e['temp']} K")
        ax.plot(e["angles"], e["vels"], 'o', color=custom_colors.get(int(e["temp"])) , label=f"{e['temp']} K")

    ax.set_title(r'Applied $\tau_{xz}$=' + s + " GPa", fontsize=12)
    ax.set_xlabel(r'Misorientation Angle [degrees]', fontsize=12)
    ax.grid(True, linestyle=':')
    ax.set_xlim(0, 180)
    ax.set_ylim(0)

# shared y-label
axes[0].set_ylabel(r'Velocity Magnitude [m/s]', fontsize=12)

# one shared legend (top)
handles, labels = axes[0].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, title="Temperature", loc='upper center',
               ncol=min(5, len(labels)), frameon=False, bbox_to_anchor=(0.5, 1.02))
fig.tight_layout(rect=(0, 0, 1, 0.95))

# Save a single figure
out_name = "Velocity_vs_Misorientation_5panels_byStress.pdf"
fig.savefig(out_name, dpi=1000, bbox_inches='tight')
print(f"Figure saved as: {out_name}")