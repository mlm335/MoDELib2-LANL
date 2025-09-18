import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------
# -- Fixed vector v and compute angle in degrees
# -------------------------------------------------------------------------
v = np.array([1, 1, -1], dtype=float)
def angle_in_degrees(u, w):
    """
    Returns the angle in DEGREES between vectors u and w.
    """
    dot_uw = np.dot(u, w)
    norm_u = np.linalg.norm(u)
    norm_w = np.linalg.norm(w)
    if norm_w == 0.0:
        return np.nan  # undefined angle for the zero vector
    cos_theta = dot_uw / (norm_u * norm_w)
    # Numerical clamp to [-1,1] to handle floating precision
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta_rad = np.arccos(cos_theta)
    return np.degrees(theta_rad)

# -------------------------------------------------------------------------
# -- 2D grid (a, b). Define w=(a,b,-a). Compute angle for each point.
# -------------------------------------------------------------------------
a_vals = np.arange(-6, 7)
b_vals = np.arange(-6, 7)
A, B = np.meshgrid(a_vals, b_vals, indexing='xy')
Theta_deg = np.zeros_like(A, dtype=float)
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        a = A[i, j]
        b = B[i, j]
        w = np.array([a, b, -a], dtype=float)
        Theta_deg[i, j] = angle_in_degrees(v, w)

# -------------------------------------------------------------------------
# -- Contour plot of angle distribution over (a, b).
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 10))
contours = ax.contourf(A, B, Theta_deg, levels=30, cmap='viridis')
cbar = plt.colorbar(contours)
cbar.set_label('Angle (degrees)')
ax.set_xlabel('a')
ax.set_ylabel('b')
ax.set_xlim(-6,6)
ax.set_ylim(-6,6)
ax.set_title(
    "Angle between v=(1,1,-1) and w=(a,b,-a)\n"
    "Overlaying lines b=a+k for multiple k."
)

# -------------------------------------------------------------------------
# -- Overlay lines b = a + k for k=1,2,3 (cross-product = k*(1,0,1)).
# -------------------------------------------------------------------------
k_values = [-3, -2, -1, 1, 2, 3]  # Choose any integer values you like
colors = ['r', 'orange', 'g', 'b', 'purple', 'k']  # One color per k, can add more if needed
sample_a_vals = [-2,-1, 0, 1,2]
for k, color in zip(k_values, colors):
    # Plot the line b = a + k
    ax.plot(a_vals, a_vals + k, color=color, linestyle='-', linewidth=2,
            label=f'Cross = ({k},0,{k}), i.e. b=a+{k}')
    # Annotate integer points along this line
    for a_pt in sample_a_vals:
        b_pt = a_pt + k
        c_pt = -a_pt
        w_pt = np.array([a_pt, b_pt, c_pt], dtype=float)
        theta_pt = angle_in_degrees(v, w_pt)
        # Mark the point with a dot
        ax.plot(a_pt, b_pt, 'o', color=color, markersize=6)
        # Annotate: "(a,b,c); angle=xx.x"
        label_text = f"({a_pt},{b_pt},{c_pt}); \n {theta_pt:.1f}°"
        # Offset the text slightly
        ax.text(a_pt-0.3, b_pt + 0.2, label_text, color=color, fontsize=8)

ax.legend(loc='best')
plt.savefig('/Users/matthewmaron/Desktop/LatticeVectos.pdf',dpi=1000)
plt.show()
