import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from matplotlib import rc

# LaTeX & font config
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['font.family'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
rc('text', usetex=True)
rc('font', family='serif')

# -------------------------------
# Load FULL trained pipeline
# -------------------------------
# loaded_model = joblib.load("SVR_model.pkl")
loaded_model = joblib.load("GBR_model.pkl")

# -------------------------------
# Utility Functions
# -------------------------------
def compute_missorientation(tangent):
    burgers = np.array([1, 1, -1]) / np.linalg.norm([1, 1, -1])
    cos_theta = np.clip(np.dot(burgers, tangent) / (np.linalg.norm(burgers) * np.linalg.norm(tangent)), -1, 1)
    return np.round(np.degrees(np.arccos(cos_theta)), 1)

def rotate_cartesian(C2G, v):
    return C2G @ v

def rotate_stress(resolved_shear, C2G):
    sigma_crystal = (resolved_shear / np.sqrt(6.0)) * np.array([[2,1,0],[1,0,1],[0,1,-2]], dtype=float)
    return C2G @ sigma_crystal @ C2G.T

def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def compute_rotation_matrix(x_dir, z_dir):
    z_target = normalize(z_dir)
    x_dir = normalize(x_dir - np.dot(x_dir, z_target) * z_target)
    y_new = normalize(np.cross(z_target, x_dir))
    return np.row_stack((x_dir, y_new, z_target))

def sort_eigenvalues_and_vectors(eigenvalues, eigenvectors):
    idx = eigenvalues.argsort()[::-1]
    return eigenvalues[idx], eigenvectors[:, idx]

# -------------------------------
# Mobility Solver Class
# -------------------------------
class MobilitySolver:
    def velocityPy(self, stress, burgers, tangent, normal, temp):
        miss_angle = compute_missorientation(tangent)
        eig_val, eig_vec = np.linalg.eig(stress)
        eig_val, eig_vec = sort_eigenvalues_and_vectors(eig_val, eig_vec)
        projection = eig_vec.T @ tangent

        # Use consistent order for columns
        features = pd.DataFrame([{
            'e1': eig_val[0],
            'e2': eig_val[1],
            'e3': eig_val[2],
            'anev1': abs(projection[0]),
            'anev2': abs(projection[1]),
            'anev3': abs(projection[2]),
            'temp': temp,
            'miss': miss_angle
        }])

        velocity_pred = loaded_model.predict(features)
        return velocity_pred[0]

# -------------------------------
# Set up Simulation
# -------------------------------
solver = MobilitySolver()

x1_lattice = np.array([1, 4, -1])   # tangent
x1_lattice = np.array([1, 1, -1])   # tangent
x3_lattice = np.array([1, 0,  1])   # normal
C2G = compute_rotation_matrix(x1_lattice, x3_lattice)

tangentV = rotate_cartesian(C2G, x1_lattice)
burgersV = rotate_cartesian(C2G, np.array([1, 1, -1]))
normalV  = rotate_cartesian(C2G, x3_lattice)

temperatures = [100, 300, 500, 700, 1000]
stress_values = np.linspace(0.0, 1.0, 1000)

results = {temp: [] for temp in temperatures}

# -------------------------------
# Run Validation Loop
# -------------------------------
for temp in temperatures:
    for stress_mag in stress_values:
        sigma = rotate_stress(stress_mag, C2G)
        velocity = solver.velocityPy(sigma, burgersV, tangentV, normalV, temp)
        results[temp].append(velocity)

missorientation = compute_missorientation(x1_lattice)

# -------------------------------
# Plot Results
# -------------------------------
cmap = plt.get_cmap("viridis")
color_map = {temp: cmap(i / len(temperatures)) for i, temp in enumerate(temperatures)}

plt.figure(figsize=(8, 8))
for temp, velocities in results.items():
    plt.plot(stress_values, velocities, color=color_map[temp], label=f'{temp}K')
plt.xlabel('Resolved Shear Stress [GPa]', fontsize=16)
plt.ylabel('Predicted Velocity [A/ps]', fontsize=16)
plt.title(f'Missorientation Angle: {missorientation}°', fontsize=20)
plt.tick_params(axis='both', which='major', direction='in', labelsize=14)
plt.tick_params(axis='both', which='minor', direction='in', labelsize=14)
plt.legend()
plt.grid()
plt.tight_layout()
# plt.savefig(f'Validation_{missorientation}degrees_SVR.pdf', dpi=1000)
plt.show()
