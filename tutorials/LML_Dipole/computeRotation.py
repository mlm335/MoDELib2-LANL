import numpy as np

def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def compute_rotation_matrix(x_dir):
    # Define fixed z-axis direction
    z_target = np.array([1, 0, 1], dtype=float)
    z_target = normalize(z_target)  # Ensure it's a unit vector
    # Make sure the input x-direction is orthogonalized with respect to z_target
    x_dir = np.array(x_dir, dtype=float)
    x_dir = x_dir - np.dot(x_dir, z_target) * z_target  # Remove z_target component
    x_dir = normalize(x_dir)  # Normalize
    # Compute the new y-axis using the cross product to ensure orthogonality
    y_new = np.cross(z_target, x_dir)
    y_new = normalize(y_new)
    # Stack the vectors as columns of the rotation matrix
    R = np.row_stack((x_dir, y_new, z_target))
    return R

# Example input vector for the x-direction
x_input = [1, -2, -1]  # Change this as needed
# Compute the rotation matrix
rotation_matrix = compute_rotation_matrix(x_input)


# Dyadic Product for the Stress Tensor Local Direction
resolved_shear = 0.1
sigma_crystal = (resolved_shear / np.sqrt(6.0)) * np.array([
    [2, 1, 0],
    [1, 0, 1],
    [0, 1,-2]
], dtype=float)
# Global Stress Tensor to keep resolved shear the same
sigma_global = rotation_matrix @ sigma_crystal @ rotation_matrix.T

# Print result
print("Rotation Matrix=\n",rotation_matrix)
print("\nSigma_crystal =\n", sigma_crystal)
print("\nSigma_global = R * sigma_crystal * R^T =\n", np.round(sigma_global,3))



