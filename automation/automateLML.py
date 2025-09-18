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

# -- Inputs to simulations
X1_Direction = [ [1,1,-1],[2,3,-2],[1,2,-1],[1,3,-1],[1,4,-1],[0,1,0],
                 [-1,2,1],[-1,1,1],[-2,1,2],[-1,0,1],[-2,1,2],[-3,-2,3] ]
Resolved_Shear = [0.0122, 0.0091463, 0.0061, 0.0031, 0.00122]
Temperature = [100, 300, 500, 700, 1000]

folder_to_run_in = 'LML_Dipole'

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


b_unit = normalize([1,1,-1])  
n_unit = normalize([1,0,1]) 
# -- Iterate through all combinations
for x1 in X1_Direction:
    x1_unit = normalize(x1)
    R = compute_rotation_matrix(x1_unit,n_unit)
    formatted_R_matrix = [[f"{val:.15f}" for val in row] for row in R]

    for shear in Resolved_Shear:
        sigma_crystal = shear* (np.outer(b_unit, n_unit) + np.outer(n_unit, b_unit))
        sigma_global = R @ sigma_crystal @ R.T
        sigma_voigt = stress_to_voigt(sigma_global)
        formatted_sigma = [f"{val:.6f}" for val in sigma_voigt]
        stress_array = np.array([float(s) for s in formatted_sigma])
        
        for temp in Temperature:
            folder_name = f"X1_{x1[0]}{x1[1]}{x1[2]}_S{shear:.4f}_T{temp}"
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
            pf.grain1globalX1=np.array(x1)     # global x1 axis. Overwritten if alignToSlipSystem0=true
            pf.grain1globalX3=np.array([1,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
            pf.boxEdges=np.array([x1,np.cross([1,0,1],x1),[1,0,1]]) # i-throw is the direction of i-th box edge
            pf.boxScaling=np.array([80,160,160])# length of box edges in Burgers vector units
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
            setInputVariable(DD_file_path,'Nsteps','5000')
            setInputVariable(DD_file_path, 'remeshFrequency', '0')
            setInputVariable(DD_file_path, 'Lmin', '20')
            setInputVariable(DD_file_path, 'Lmax', '25')
            setInputVariable(DD_file_path, 'alphaLineTension', '1')
            setInputVariable(DD_file_path, 'timeSteppingMethod', 'fixed')  # fixed or adaptive
            setInputVariable(DD_file_path, 'dtMax', '8')
            setInputVariable(DD_file_path, 'dxMax', '1.4')
            setInputVector(DD_file_path, 'periodicImageSize', np.array([8, 8, 8]),'' )
            setInputVariable(DD_file_path, 'EwaldLengthFactor', '4')
            setInputVariable(DD_file_path, 'quadPerLength', '0.001')


            # -- Microstructure File
            microstructureFile1='periodicDipoleIndividual.txt';
            microstructureFileTemplate='../Library/Microstructures/'+microstructureFile1;
            print("\033[1;32mCreating  microstructureFile\033[0m")
            microstructure_file_path1 = os.path.join(input_directory, microstructureFile1)
            shutil.copy2(microstructureFileTemplate, microstructure_file_path1) # target filename is /dst/dir/file.ext
            setInputVariable(microstructure_file_path1, "slipSystemIDs", "0")
            setInputVariable(microstructure_file_path1, "exitFaceIDs", "3")
            setInputVariable(microstructure_file_path1, "nodesPerLine", "1")
            setInputVector(microstructure_file_path1, "dipoleCenters", np.array([0, 0, 0]),'')
            setInputVariable(microstructure_file_path1, "dipoleHeights", "100")
            setInputVariable(microstructure_file_path1, "glideSteps", "40")


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
