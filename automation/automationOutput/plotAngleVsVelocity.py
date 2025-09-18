import math, sys, string, os
from matplotlib import rc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
print(os.environ['PATH'])
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['font.family'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
rc('text', usetex=True)
rc('font', family='serif')
sys.path.append('../../lib')

# -- Load Data
csv_filename = "simulation_results.csv"
df = pd.read_csv(csv_filename)
unique_stresses = np.unique(df["Stress [GPa]"])

# -- Colormap for temperatures
cmap = plt.get_cmap("viridis")
unique_temperatures = np.unique(df["Temperature [K]"])
color_map = {temp: cmap(i / len(unique_temperatures)) for i, temp in enumerate(unique_temperatures)}

# Plots for each stress value
for stress in unique_stresses:
    subset = df[df["Stress [GPa]"] == stress]
    fig, ax = plt.subplots(figsize=(8, 6))
    for temp in unique_temperatures:
        temp_subset = subset[subset["Temperature [K]"] == temp]
        # -- Sort by misorientation angle
        temp_subset = temp_subset.sort_values(by="Misorientation Angle [°]")
        if not temp_subset.empty:
            ax.plot(temp_subset["Misorientation Angle [°]"],
                    temp_subset["Time Average Velocity [Ang/ps]"],
                    label=f"{temp} K", color=color_map[temp],
                    linestyle="--", marker='o', markersize=6)

    ax.set_xlabel("Misorientation Angle [°]",fontsize=16)
    ax.set_ylabel("Time Average Velocity [Ang/ps]",fontsize=16)
    ax.set_title(f"Resolved Shear Stress = {stress} GPa",fontsize=16)
    ax.legend(title="Temperature [K]", fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)

    # Save each plot with stress in filename
    plot_filename = f"Velocity_vs_Misorientation_Stress_{stress:.2f}GPa.pdf"
    plt.savefig(plot_filename, dpi=1000)
    print(f"Saved {plot_filename}")

plt.show()
