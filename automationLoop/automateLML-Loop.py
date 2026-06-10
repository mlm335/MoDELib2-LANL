import subprocess, os, sys, pathlib, json, time
import numpy as np
sys.path.append('lib')
from write_data import *
from read_data import *
from generateInputFiles import *
from fractions import Fraction
import shutil
sys.path.append("../python/")
from modlibUtils import *


# -- Inputs to simulations
Resolved_Shear = [0.0122, 0.0091463, 0.0061, 0.0031, 0.00122]
Temperature = [100, 300, 500, 700, 1000]

folder_to_run_in = 'LML_Loop'

# -- Working directories
current_dir = str(pathlib.Path().resolve())
output_directory = os.path.join(current_dir, "automationOutput")
os.makedirs(output_directory, exist_ok=True)
folder_directory = os.path.join(current_dir, "../tutorials/", folder_to_run_in)
input_directory = os.path.join(folder_directory, "inputFiles")
tools_directory = os.path.join(current_dir, "../build/tools/")
MicrostructureGenerator = os.path.join(tools_directory, "MicrostructureGenerator/")
DDomp = os.path.join(tools_directory, "DDomp/")
original_dir = os.getcwd()


for shear in Resolved_Shear:
    stress_array = np.array([0.0 ,0.0, 0.0, 0.0, 0.0, shear])  # only sigma_xz is non-zero
    
    for temp in Temperature:
        folder_name = f"S{shear:.4f}_T{temp}"
        output_path = os.path.join(output_directory, folder_name)
        os.makedirs(output_path, exist_ok=True)
        print(f"Running simulation for {folder_name}")

        os.chdir(original_dir)  # <- restore path before relative-copy operations for next iteration

        # -- Material File
        materialFile='Fe_320.txt';
        materialFileTemplate='../Library/Materials/'+materialFile;
        print("\033[1;32mCreating  materialFile\033[0m")
        shutil.copy2(materialFileTemplate, os.path.join(input_directory, materialFile))

        # -- Polycrystal File
        meshFile = 'unitCube24.msh'
        meshFileTemplate = '../Library/Meshes/' + meshFile
        shutil.copy2(meshFileTemplate, os.path.join(input_directory, meshFile))
        os.chdir(folder_directory) # <- relative-copy operations 
        pf=PolyCrystalFile(materialFile);
        pf.absoluteTemperature=temp;
        pf.meshFile=meshFile
        pf.grain1globalX1=np.array([1,1,-1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
        pf.grain1globalX3=np.array([1,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
        pf.boxEdges=np.array([[1,1,-1],np.cross([1,0,1],[1,1,-1]),[1,0,1]]) # i-throw is the direction of i-th box edge
        pf.boxScaling=np.array([200,200,100])# length of box edges in Burgers vector units
        pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
        pf.periodicFaceIDs=np.array([-1])
        pf.write('inputFiles')

        os.chdir(original_dir)  # <- restore path after relative-copy operations 


        # -- Elastic Deformation File
        elasticDeformatinoFile='ElasticDeformation.txt';
        elasticDeformatinoFileTemplate='../Library/ElasticDeformation/'+elasticDeformatinoFile;
        print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
        elastic_file_path = os.path.join(input_directory, "ElasticDeformation.txt")
        shutil.copy2(elasticDeformatinoFileTemplate,elastic_file_path)
        setInputVector(elastic_file_path, "ExternalStress0", stress_array,'')


        # -- DD Parameters File 
        DDfile='DD.txt'
        DDfileTemplate='../Library/DislocationDynamics/'+DDfile
        print("\033[1;32mCreating  DDfile\033[0m")
        DD_file_path = os.path.join(input_directory, DDfile)
        shutil.copy2(DDfileTemplate, DD_file_path)
        setInputVariable(DD_file_path,'useFEM','0')
        setInputVariable(DD_file_path,'useDislocations','1')
        setInputVariable(DD_file_path,'useInclusions','0')
        setInputVariable(DD_file_path,'useElasticDeformation','1')
        setInputVariable(DD_file_path,'useClusterDynamics','0')
        setInputVariable(DD_file_path,'Nsteps','300')
        setInputVariable(DD_file_path, 'remeshFrequency', '1')
        setInputVariable(DD_file_path, 'Lmin', '30')
        setInputVariable(DD_file_path, 'Lmax', '50')
        setInputVariable(DD_file_path, 'alphaLineTension', '0')
        setInputVariable(DD_file_path,'computeDDinteractions','0')  # DD interactions
        setInputVariable(DD_file_path, 'timeSteppingMethod', 'fixed')  # fixed or adaptive
        setInputVariable(DD_file_path, 'dtMax', '200')
        setInputVariable(DD_file_path, 'dxMax', '1.0')
        setInputVariable(DD_file_path,'use_velocityFilter','0') # don't filter velocity if noise is enabled
        setInputVariable(DD_file_path,'use_stochasticForce','0') # Langevin thermal noise enabled
        setInputVector(DD_file_path, 'periodicImageSize', np.array([8, 8, 8]),'' )
        setInputVariable(DD_file_path, 'EwaldLengthFactor', '4')
        setInputVariable(DD_file_path, 'quadPerLength', '0.0001')
        setInputVariable(DD_file_path,'outputFrequency','1')  # output frequency
        setInputVariable(DD_file_path,'outputQuadraturePoints','1')  # output quadrature data
        setInputVariable(DD_file_path,'glideSolverType','Galerkin')  # type of glide solver, or none
        setInputVariable(DD_file_path,'climbSolverType','none')  # type of clim solver, or none

        # -- Microstructure File        
        microstructureFile1='shearLoopsIndividual.txt';
        microstructureFileTemplate='../Library/Microstructures/'+microstructureFile1;
        print("\033[1;32mCreating  microstructureFile\033[0m")
        microstructure_file_path1 = os.path.join(input_directory, microstructureFile1)
        shutil.copy2(microstructureFileTemplate, microstructure_file_path1) # target filename is /dst/dir/file.ext
        setInputVector(microstructure_file_path1,'slipSystemIDs',np.array([1]),'')
        setInputVector(microstructure_file_path1,'loopRadii_SI',np.array([10e-9]),'')
        setInputVector(microstructure_file_path1,'loopSides',np.array([28]),'')
        setInputVector(microstructure_file_path1,'loopCenters',np.array([0, 0, 0]),'')


        print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
        initial_micro_path = os.path.join(input_directory, "initialMicrostructure.txt")
        with open(initial_micro_path, "w") as f:
                f.write(f'microstructureFile={microstructureFile1};\n')


        # -- Run the simulation
        time.sleep(1)
        os.chdir(folder_directory)
        os.system('/Users/matthewmaron/bashScripts/allrmEvlF.sh')
        time.sleep(1)
        os.chdir(MicrostructureGenerator)
        os.system(f'./microstructureGenerator {folder_directory}')
        time.sleep(1)
        os.chdir(DDomp)
        os.system(f'./DDomp {folder_directory}')
        time.sleep(1)
        os.chdir(folder_directory)
        

        # -- Save output
        evl_path = generate_folder(output_path, "evl")
        F_path = generate_folder(output_path, "F")
        input_path = generate_folder(output_path, "inputFiles")
        copy_folder_contents('F', F_path)
        copy_folder_contents('evl', evl_path)
        copy_folder_contents('inputFiles', input_path)
        print(f"Simulation for {folder_name} completed and saved.")
