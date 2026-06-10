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

# -- Loop and Quadrature Point Data
t, b = getLoopData(dataFolder)
t2, tangent, velocity = getQPData(dataFolder)

# -- Misorientation Angle
burgers = np.array(b[0])
miss_angle = []
for k in range(len(tangent[0])):
    QPtangent = np.array(tangent[0][k])
    miss_product = np.dot(burgers, QPtangent)
    magnitude_a = np.linalg.norm(burgers)
    magnitude_b = np.linalg.norm(QPtangent)
    cos_theta = miss_product / (magnitude_a * magnitude_b)
    angle_radians = np.arccos(cos_theta)
    miss_angle.append(np.degrees(angle_radians))

# -- Velocity Magnitude
VelocityMagnitudes = np.linalg.norm(velocity[0], axis=1)

# -- Plot for this folder
plt.plot(
    miss_angle,
    VelocityMagnitudes * cs,
    'o',
    label=f"{int(Temperature)} K"
)

# Beautification
plt.legend()
plt.title(r'Applied $\tau_{xz}$=' + str(np.round(Stress0, 2)) + " GPa", fontsize=12)
plt.xlabel(r'Misorientation Angle [degrees]',fontsize=12)
plt.ylabel(r'Velocity Magnitude [m/s]',fontsize=12)
plt.xlim(0, 180)
plt.ylim(0)
plt.grid(True, linestyle=':')

# Save Figure
stress_str = f"{np.round(Stress0, 2)}"
filename = f"Velocity_vs_Misorientation_Stress_{stress_str}_GPa.pdf"
plt.savefig(filename, dpi=1000, bbox_inches='tight')
print(f"Figure saved as: {filename}")

# Show Figure
#plt.show()
