import sys
sys.path.append("../../python/")
from modlibUtils import *

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
setInputVariable('inputFiles/'+DDfile,'useFEM','0')
setInputVariable('inputFiles/'+DDfile,'useDislocations','1')
setInputVariable('inputFiles/'+DDfile,'useInclusions','0')
setInputVariable('inputFiles/'+DDfile,'useElasticDeformation','1')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamics','0')

setInputVariable('inputFiles/'+DDfile,'Nsteps','100')  # number of simulation steps
setInputVector('inputFiles/'  +DDfile,'periodicImageSize',np.array([8,8,8]),'')
setInputVariable('inputFiles/'+DDfile,'EwaldLengthFactor','4')
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','fixed') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','2')
setInputVariable('inputFiles/'+DDfile,'dxMax','1') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0') # don't filter velocity 
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'useSubCycling','0')

setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'quadPerLength','0.1')
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','0')  # DD interactions
setInputVariable('inputFiles/'+DDfile,'Lmin','10')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','15')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','10')
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'computeDDinteractions','1')  # DD interactions


# Make a local copy of material file, and modify that copy if necessary
materialFile='W.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'dislocationMobilityType','default')
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','full<111>{110}')


# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.01]),'applied stress')


# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile

# pf.grain1globalX1=np.array([1,1,-1])  # Pure EDGE
# pf.grain1globalX1=np.array([1,-2,-1])  # Pure SCREW
# pf.grain1globalX1=np.array([1,0,-1])    # Mixed Theta about 35 degrees
# pf.grain1globalX1=np.array([2,3,-2])  # 11.4 degrees bxt

x1=[1,1,-1]
x3=[1,0,1]
x2=np.cross(x3,x1)

pf.grain1globalX1=np.array(x1)
pf.grain1globalX3=np.array(x3)    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([x1,x2,x3]) # i-throw is the direction of i-th box edge

# basisLength = np.linalg.norm([1,1,-1])
# f1 = 400/2*basisLength / np.linalg.norm(x1)
# f2 = 100*basisLength / np.linalg.norm(x2)
# f3 = 200/2 * basisLength / np.linalg.norm(x3)

pf.boxScaling=np.array([80,160,160]) # must be a vector of integers
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')


# make a local copy of microstructure file, and modify that copy if necessary
microstructureFile1='periodicDipoleIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVariable('inputFiles/'+microstructureFile1,'slipSystemIDs','0')
setInputVariable('inputFiles/'+microstructureFile1,'exitFaceIDs','3')
setInputVector('inputFiles/'+microstructureFile1,'dipoleCenters',np.array([0.0,0.0,0.0]),'')
setInputVariable('inputFiles/'+microstructureFile1,'nodesPerLine','20')
setInputVariable('inputFiles/'+microstructureFile1,'dipoleHeights','100')
setInputVariable('inputFiles/'+microstructureFile1,'glideSteps','40.0')


print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
#    initialMicrostructureFile.write('microstructureFile='+microstructureFile2+';\n')
