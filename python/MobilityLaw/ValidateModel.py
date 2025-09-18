import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from matplotlib import rc
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['font.family'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
rc('text', usetex=True)
rc('font', family='serif')

# Load the model
# loaded_model = joblib.load('/Users/matthewmaron/Documents/MoDELib2-new/python/MobilityLaw/GPR_model.pkl')
# loaded_model = joblib.load('/Users/matthewmaron/Documents/MoDELib2-new/python/MobilityLaw/SVR_model.pkl')
# loaded_model = joblib.load('/Users/matthewmaron/Documents/MoDELib2-new/python/MobilityLaw/KRR_model.pkl')

# loaded_model = joblib.load('/Users/matthewmaron/Documents/MoDELib2-new/python/MobilityLaw/GBR_model.pkl')
loaded_model = joblib.load('/Users/matthewmaron/Documents/MoDELib2-new/python/MobilityLaw/MLP_model.pkl')

# ----------------------------------------------------------------- #
def compute_missorientation(tangent):
    # -- Misorientation Angle
    burgers = np.array([1,1,-1]) / np.linalg.norm(np.array([1,1,-1]))
    miss_product = np.dot(burgers, tangent)
    magnitude_a = np.linalg.norm(burgers)
    magnitude_b = np.linalg.norm(tangent)
    cos_theta = np.clip(miss_product / (magnitude_a * magnitude_b), -1, 1)
    angle_radians = np.arccos(cos_theta)
    MissorientationAngle = np.round(np.degrees(angle_radians),1)
    return MissorientationAngle
    
def rotate_cartesian(C2G,v):
    return C2G @ v
    
def rotate_stress(resolved_shear, C2G):
    sigma_crystal = (resolved_shear / np.sqrt(6.0)) * np.array([
        [2, 1, 0],
        [1, 0, 1],
        [0, 1,-2]
    ], dtype=float)
    sigma_global = C2G @ sigma_crystal @ C2G.T
    return sigma_global
        
def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def compute_rotation_matrix(x_dir,z_dir):
    z_target = np.array(z_dir, dtype=float)
    z_target = normalize(z_target)
    x_dir = np.array(x_dir, dtype=float)
    x_dir = x_dir - np.dot(x_dir, z_target) * z_target
    x_dir = normalize(x_dir)
    y_new = np.cross(z_target, x_dir)
    y_new = normalize(y_new)
    return np.row_stack((x_dir, y_new, z_target))
    
def sort_eigenvalues_and_vectors(eigenvalues, eigenvectors):
    idx = eigenvalues.argsort()[::-1]  # Sort in descending order
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues, eigenvectors

class MobilitySolver():
    def velocityPy(self, stress, burgers, tangent, normal, temp):
    
        miss_product = np.dot(burgers, tangent)
        magnitude_a = np.linalg.norm(burgers)
        magnitude_b = np.linalg.norm(tangent)
        cos_theta = np.clip(miss_product / (magnitude_a * magnitude_b), -1, 1)
        angle_radians = np.arccos(cos_theta)
        miss_angle = np.degrees(angle_radians)
        
        eigValue, eigVector = np.linalg.eig(stress)
        eigValue, eigVector = sort_eigenvalues_and_vectors(eigValue, eigVector)
        projection = eigVector.T @ tangent
        
        features = pd.DataFrame({
            'e1': eigValue[0],
            'e2': eigValue[1],
            'e3': eigValue[2],
            # 'anev1': abs(projection[0]),
            # 'anev2': abs(projection[1]),
            # 'anev3': abs(projection[2]),
            'temp': temp,
            'miss': miss_angle
        }, index=[0])

        # print(features)

        velocity_pred = loaded_model.predict(features)
    
        # print(f"Predicted Velocity: {velocity_pred[0]} A/ps")

        return velocity_pred[0]
# ----------------------------------------------------------------- #

# -- Instance of MobilitySolver
solver = MobilitySolver()

# -- Lattice Directions and C2G1
# x1_lattice = np.array([1,4,-1]) # tangent for mixed
x1_lattice = np.array([-1,2,1]) # tangent for edge
# x1_lattice = np.array([1,1,-1]) # tangent for screw

x3_lattice = np.array([1,0,1]) # normal
C2G = compute_rotation_matrix(x1_lattice,x3_lattice)

# -- Cartesian
tangentV = rotate_cartesian(C2G,x1_lattice)
burgersV = rotate_cartesian(C2G,np.array([1,1,-1]))
normalV = rotate_cartesian(C2G,x3_lattice)

temperatures = [100, 300, 500, 700, 1000]  # K
stress_values = np.linspace(0.01, 1, 1000) # GPa


# -- Store results
results = {temp: [] for temp in temperatures}

# -- Loop
for temp in temperatures:
    for stress_mag in stress_values:
        Stress0 = rotate_stress(stress_mag,C2G)
 
        velocity = solver.velocityPy(Stress0, burgersV, tangentV, normalV, temp)
        results[temp].append(velocity)

# -- Missorientation
missorientation = compute_missorientation(x1_lattice)

# -- Plot
# -- Colormap for temperatures
cmap = plt.get_cmap("viridis")
color_map = {temp: cmap(i / len(temperatures)) for i, temp in enumerate(temperatures)}
plt.figure(figsize=(8, 8))
for temp, velocities in results.items():
    plt.plot(stress_values, velocities, color=color_map[temp], label=f'{temp}K')
plt.xlabel('Resolved Shear Stress [GPa]',fontsize=16)
plt.ylabel('Predicted Velocity [A/ps]',fontsize=16)
plt.title(f'Missorientation Angle: {missorientation}°',fontsize=20)
plt.tick_params(axis = 'both', which = 'major', direction='in', bottom=True, left=True, top=True, right=True, labelsize = 14)
plt.tick_params(axis = 'both', which = 'minor', direction='in', bottom=True, left=True, top=True, right=True, labelsize = 14)
plt.legend()
plt.grid()
plt.savefig(f'Validation_{missorientation}degrees_MLP.pdf',dpi=1000)
plt.show()

