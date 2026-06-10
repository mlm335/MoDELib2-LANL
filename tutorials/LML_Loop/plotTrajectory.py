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
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['font.family'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
rc('text', usetex=True)
rc('font', family='serif')
sys.path.append('../../python/lib')
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
    pos=[];
    velo=[];
    n=0;
    for runID in F[:,0]:
        evl=readEVLtxt(folderName+'/evl/evl_'+str(int(runID)))
        b.append(evl.loopsBurger);
        pos.append(evl.nodesPos)
        velo.append(evl.nodesV);
        n=n+1;
    return F[:,1], b, pos, velo;

###########################################################################
###########################################################################
###########################################################################


dataFolder = str(pathlib.Path().resolve())

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

t, loopBurger, nodePos, nodeVelo = getLoopData(dataFolder)


def plot_loops_trajectory(t, nodePos, n_points=50, cmap_name="viridis"):

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    n_total = len(t)
    if n_points is None or n_points >= n_total:
        idxs = np.arange(n_total)
    else:
        # pick evenly spaced indices
        idxs = np.linspace(0, n_total-1, n_points, dtype=int)

    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0, 1, len(idxs)))

    for color, k in zip(colors, idxs):
        nodes = np.array(nodePos[k])
        xy = nodes[:, :2]
        ax.scatter(xy[:,0], xy[:,1], np.full(xy.shape[0], t[k])*t_dd2SI*1e12,
                   s=10, color=color, label=f"t={t[k]:.2f}")

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
    plt.show()

# ---- usage ----
plot_loops_trajectory(t[0:100], nodePos[0:100], n_points=3)   # only 10 time step
