# Some useful functions
import numpy as np
import math, sys, string, os
from matplotlib import pyplot as plt
from matplotlib import rc
import pandas as pd
import matplotlib.patches as patches
from os import listdir
from os.path import isfile, join
from mpl_toolkits.mplot3d import Axes3D  # needed for 3D
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
    nID=[];
    b=[];
    pos=[];
    velo=[];
    n=0;
    for runID in F[:,0]:
        evl=readEVLtxt(folderName+'/evl/evl_'+str(int(runID)))
        nID.append(evl.nodesID);
        b.append(evl.loopsBurger);
        pos.append(evl.nodesPos)
        velo.append(evl.nodesV);
        n=n+1;
    return F[:,1], nID, b, pos, velo;


def plot_loops_trajectory(directory, folderinfo, t, nodeID, nodePos, n_points=50, cmap_name="viridis"):

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    n_total = len(t)
    if n_points is None or n_points >= n_total:
        idxs = np.arange(n_total)
    else:
        idxs = np.linspace(0, n_total-1, n_points, dtype=int)  # pick evenly spaced indices

    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0, 1, len(idxs)))

    # ---- Plot the Nodes ----
    for color, k in zip(colors, idxs):
        nodes = np.array(nodePos[k])
        xy = nodes[:, :2]
        z_val = float(t[k]) * t_dd2SI * 1e12      # ps (uses your global t_dd2SI)
        ax.scatter(xy[:,0], xy[:,1], np.full(xy.shape[0], z_val),
                s=10, color=color, label=f"t={t[k]:.2f}")

        # ---- connect nodes by ID order: 0-1, 1-2, ..., (N-1)-0 ----
        N = xy.shape[0]
        if N >= 2:
            for i in range(N):
                j = (i + 1) % N  # wrap to close the loop
                ax.plot( [xy[i,0], xy[j,0]], [xy[i,1], xy[j,1]], [z_val,   z_val], '-', linewidth=1.0, color=color)

    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    # ---- Hide X and Y axes completely ----
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.xaxis.set_pane_color((1,1,1,0))   # transparent background
    ax.yaxis.set_pane_color((1,1,1,0))
    ax.zaxis.set_pane_color((1,1,1,0))
    ax.xaxis.line.set_color((1,1,1,0))   # hide axis line
    ax.yaxis.line.set_color((1,1,1,0))
    ax.zaxis.line.set_color((1,1,1,0))
    ax.xaxis._axinfo['grid']['color'] = (1,1,1,0)  # hide grid
    ax.yaxis._axinfo['grid']['color'] = (1,1,1,0)
    ax.zaxis._axinfo['grid']['color'] = (1,1,1,0)

    ax.set_zlabel("time [ps]")
    ax.view_init(elev=25, azim=-60)
    plt.tight_layout()
    T = list(folderinfo)[0]
    Stress0 = list(folderinfo)[1]
    plt.savefig(f"{directory}/TrajectoriesPlots/{T}_{Stress0}.pdf", dpi=1000)
    # plt.show()
    
###########################################################################
###########################################################################
###########################################################################
current_dir = str(pathlib.Path().resolve())
parent_directory = str(pathlib.Path().resolve())+'/automationOutput'
print("Parent Directory = ", parent_directory)
dataFolders = [f for f in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, f))]

for folder in dataFolders:
    dataFolder = os.path.join(parent_directory, folder)


    # -- Material properties
    mat_filepath = dataFolder + '/inputFiles/Fe_320.txt'
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

    t, nodeID, loopBurger, nodePos, nodeVelo = getLoopData(dataFolder)

    folderinfo = {Temperature, np.round(Stress0,2)} # for filename
    # ---- usage ----
    plot_loops_trajectory(current_dir, folderinfo, t, nodeID, nodePos, n_points=5) 
