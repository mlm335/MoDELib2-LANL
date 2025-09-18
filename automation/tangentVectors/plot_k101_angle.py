import numpy as np
import matplotlib.pyplot as plt

def compute_angle_deg(a, b):
    """Angle (deg) between b=(1,1,-1) and t=(a,b,-a)."""
    dot = 2.0 * a + b
    norm_v = np.sqrt(3.0)
    norm_w = np.sqrt(2.0 * a * a + b * b)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_th = dot / (norm_v * norm_w)
    cos_th = np.clip(cos_th, -1.0, 1.0)
    theta = np.degrees(np.arccos(cos_th))
    return np.ma.masked_where(norm_w == 0.0, theta)


# Config
a_min, a_max = -6.0, 6.0
b_min, b_max = -6.0, 6.0
num_pts = 500

k_values = [-3, -2, -1, 0, 1, 2, 3]  # lines b = a + k to overlay

point_marker_size = 30
point_alpha = 0.95
label_fontsize = 10
label_y_offset = 0.14
show_angle_on_label = True

outfile = "angle_k101_plot.pdf"


# Build angle field
A = np.linspace(a_min, a_max, num_pts)
B = np.linspace(b_min, b_max, num_pts)
a, b = np.meshgrid(A, B)
theta_deg = compute_angle_deg(a, b)


# Plot field + lines
fig, ax = plt.subplots(figsize=(12, 10))
cont = ax.contourf(a, b, theta_deg, levels=180, cmap="viridis")
cs = ax.contour(a, b, theta_deg, levels=np.arange(0, 181, 15),
                colors="k", alpha=0.25, linewidths=0.6)
ax.clabel(cs, fmt="%d", inline=True, fontsize=7)
# Draw lines b = a + k and record colors
k_to_color = {}
for k in k_values:
    (line,) = ax.plot(A, A + k, linewidth=2.4, label=f"k={k}")
    k_to_color[k] = line.get_color()


# Plot integer points ONLY on lines, not on edges
# interior integer ranges (exclude edges exactly at min/max)
a_ints = range(int(np.ceil(a_min)) + 1, int(np.floor(a_max)))
b_ints = range(int(np.ceil(b_min)) + 1, int(np.floor(b_max)))

for ai in a_ints:
    for bi in b_ints:
        k_here = bi - ai
        if k_here not in k_values:
            continue  # only keep points that lie on one of the drawn lines
        ang = compute_angle_deg(float(ai), float(bi))
        if np.ma.is_masked(ang):
            continue
        ang_val = float(ang)
        color = k_to_color[k_here]
        ax.scatter([ai], [bi], s=point_marker_size, color=color, alpha=point_alpha, zorder=3)
        label = f"({ai},{bi},{-ai})"
        if show_angle_on_label:
            label += f"\n{ang_val:.1f}°"
        ax.text(ai, bi + label_y_offset, label,
                fontsize=label_fontsize, ha="center", va="bottom",
                color=color, zorder=3)

# Cosmetics
ax.legend(loc="upper left", framealpha=0.9, fontsize=12, title="b = a + k")
ax.set_aspect("equal", adjustable="box")
ax.set_xlim(a_min, a_max)
ax.set_ylim(b_min, b_max)
ax.set_xlabel("a",fontsize=20)
ax.set_ylabel("b",fontsize=20)
ax.tick_params(axis="both", which="major", labelsize=14)  # x and y ticks
plt.rcParams.update({
    "font.size": 20,       # base font size
    "axes.labelsize": 20,  # axis labels
    "xtick.labelsize": 16, # x tick labels
    "ytick.labelsize": 16, # y tick labels
    "legend.fontsize": 16, # legend
    "axes.titlesize": 18,  # title
})

cbar = fig.colorbar(cont, ax=ax)
cbar.set_label("Angle (degrees)")

plt.tight_layout()
plt.savefig(outfile, dpi=1000, bbox_inches="tight")
print(f"Saved figure to {outfile}")
plt.show()
