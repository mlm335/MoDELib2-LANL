import sys
sys.path.append("../../python/")
from modlibUtils import *
import numpy as np

# Create folder structure
folders=['evl','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate='../../Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
# General Settings
setInputVariable('inputFiles/'+DDfile,'useFEM','0')
setInputVariable('inputFiles/'+DDfile,'useElasticDeformation','1')
setInputVariable('inputFiles/'+DDfile,'useElasticDeformationFEM','1')
setInputVariable('inputFiles/'+DDfile,'useDislocations','1')
setInputVariable('inputFiles/'+DDfile,'useInclusions','0')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamics','0')
setInputVariable('inputFiles/'+DDfile,'Nsteps','100000')  # number of simulation steps
setInputVector('inputFiles/'+DDfile,'periodicImageSize',np.array([1, 1, 1]),'number of images for each periodic shift vector') 
setInputVariable('inputFiles/'+DDfile,'EwaldLengthFactor','1') 
# Time Stepping Settings
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','fixed') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','0.5')
setInputVariable('inputFiles/'+DDfile,'dxMax','1.234') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'velocityReductionFactor','0.75') 
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'useSubCycling','0') 
# Glide and Cllimb Solver Settings
setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')  # type of clim solver, or none
# Dislcocation Elastic Fields and Remeshing Settings
setInputVariable('inputFiles/'+DDfile,'quadPerLength','0.051')  # number of quadrature points per unit length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'coreSize','2.0')  # number of quadrature points per unit length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','6.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','100')
setInputVariable('inputFiles/'+DDfile,'Lmin','2')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','8')  # max segment length (in Burgers vector units)
# Cross Slip Junctions and Output Settings
setInputVariable('inputFiles/'+DDfile,'maxJunctionIterations','1')  
setInputVariable('inputFiles/'+DDfile,'crossSlipModel','0')  # crossSlipModel
setInputVariable('inputFiles/'+DDfile,'crossSlipAngle_deg','5')  
setInputVariable('inputFiles/'+DDfile,'outputFrequency','100')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','0')  # output quadrature data


# Make a local copy of material file, and modify that copy if necessary
materialFile='Fe_320.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','full<111>{110}')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# MD Simulation Parameters 
XLength = 20.5e-9/b_SI   # meters/b
YLength = 22.25e-9/b_SI     # meters/b
ZLength = 110.34e-9/b_SI    # meters/b
print("XLength (in b): ", XLength, " YLength (in b): ", YLength, " ZLength (in b): ", ZLength)

t_dd2SI = 1.083093316853979e-13
StrainRate = 2.5e7*t_dd2SI # in 1/s * s

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStrainRate',np.array([0.0,0.0,-StrainRate,0.0,0.0,0.0]),'applied strain rate')
setInputVector('inputFiles/'+elasticDeformatinoFile,'stiffnessRatio',np.array([0.0,0.0,1e20,0.0,0.0,0.0]),'machine stiffness')



# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
pf.grain1globalX1=np.array([3,1,0])
pf.grain1globalX3=np.array([0,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([[3,1,0],np.cross([0,0,1],[3,1,0]),[0,0,1]]) # i-throw is the direction of i-th box edge
pf.boxScaling=np.array([XLength,YLength,ZLength])# length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')

# make a local copy of microstructure file, and modify that copy if necessary
microstructureFile1='prismaticLoopsIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVector('inputFiles/'+microstructureFile1,'slipSystemIDs',np.array([5,6,11]),'slip system IDs for each loop')
setInputVector('inputFiles/'+microstructureFile1,'loopRadii_SI',np.array([1.9e-9,1.9e-9,1.9e-9]),'radii of loops in SI')
setInputVector('inputFiles/'+microstructureFile1,'glideSteps',np.array([0,0,0]),'[b] distance from sessile loop (not important, used for debug)')
setInputMatrix('inputFiles/'+microstructureFile1,'loopCenters',np.array([[0,0,0],[-XLength/4,-YLength/4,-ZLength/4],[XLength/4,YLength/4,ZLength/4]]))


print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
