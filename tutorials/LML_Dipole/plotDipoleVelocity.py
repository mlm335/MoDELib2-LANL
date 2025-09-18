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
sys.path.append('../../Results/lib')
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
    position=[]
    n=0;
    for runID in F[:,0]:
        aux=readAUXtxt(folderName+'/evl/ddAux_'+str(int(runID)))
        tangent.append(aux.tangent);
        velocity.append(aux.velocity);
        position.append(aux.position)
        n=n+1;
    return F[:,1], tangent, velocity, position;
###########################################################################
###########################################################################
###########################################################################


# Define the parent directory where your folders are located
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
C2G = get_matrix(poly_filepath, 'C2G1')

# -- Loop and Quadrature Point Data
t, b = getLoopData(dataFolder)
t2, tangent, velocity, position = getQPData(dataFolder)

# -- ElasticDeformation
NormalDirection = np.array([1,0,1])
Elastic_filepath = dataFolder + '/inputFiles/ElasticDeformation.txt'
ExternalStress0 = get_vector(Elastic_filepath, 'ExternalStress0')
StressTensor = np.array([
    [ExternalStress0[0], ExternalStress0[3], ExternalStress0[5]],  # [σ11, σ12, σ13]
    [ExternalStress0[3], ExternalStress0[1], ExternalStress0[4]],  # [σ12, σ22, σ23]
    [ExternalStress0[5], ExternalStress0[4], ExternalStress0[2]]   # [σ13, σ23, σ33]
])
Stress0 = np.dot( np.dot(StressTensor,b[0][0]),NormalDirection ) *mu0_SI*1e-9

# -- Misorientation Angle
burgers = np.array([1,1,-1])
QPtangent = np.array(C2G[0])
miss_product = np.dot(burgers, QPtangent)
magnitude_a = np.linalg.norm(burgers)
magnitude_b = np.linalg.norm(QPtangent)
cos_theta = np.clip(miss_product / (magnitude_a * magnitude_b), -1, 1)
angle_radians = np.arccos(cos_theta)
MissorientationAngle =np.degrees(angle_radians)
print("Missorientation Angle [degrees]: ", MissorientationAngle)

# -- Velocity Magnitude
average_position = [np.mean(np.linalg.norm(position[k][0:7], axis=1)) for k in range(len(velocity))]
average_velocity = [np.mean(np.linalg.norm(velocity[k], axis=1))*v_dd2SI/100 for k in range(len(velocity))]


fig1, ax1 = plt.subplots(figsize=(8, 6))
fig2, ax2 = plt.subplots(figsize=(8, 6))

F,Flabels=readFfile(dataFolder)
runID = getFarray(F,Flabels,'runID');

TimeAverageVelocity = np.mean(average_velocity)
print("Average Velocity [Ang/ps]: ", TimeAverageVelocity)

# Plot average velocity over time
ax1.plot(t2*t_dd2SI*10e12, average_velocity, label=f"{float(MissorientationAngle),float(TimeAverageVelocity) }", marker='o')
#ax1.plot(runID[0:100], average_velocity[0:100], label=f"{float(MissorientationAngle),float(TimeAverageVelocity) }", marker='o')
ax2.scatter(average_position, average_velocity, label=f"{float(MissorientationAngle),float(TimeAverageVelocity) }", marker='o')


ax1.set_xlabel("Time [ps]")
ax1.set_ylabel("Velocity [Ang/ps]")
ax2.set_xlabel("Position [b]")
ax2.set_ylabel("Velocity [Ang/ps]")
ax1.legend()
ax2.legend()

plt.show()
