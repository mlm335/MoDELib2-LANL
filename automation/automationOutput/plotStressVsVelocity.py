import math, sys, string, os
from matplotlib import rc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Font and LaTeX settings
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['font.family'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
rc('text', usetex=True)
rc('font', family='serif')
sys.path.append('../../lib')

# Load data
csv_filename = "simulation_results.csv"
df = pd.read_csv(csv_filename)

# Unique values
unique_stresses = np.unique(df["Stress [GPa]"])
unique_temperatures = np.unique(df["Temperature [K]"])
unique_misorientations = np.unique(df["Misorientation Angle [°]"])

# Pick misorientations to plot
selected_angles = sorted(unique_misorientations)[:]  

# Colormap
cmap = plt.get_cmap("viridis")
color_map = {temp: cmap(i / len(unique_temperatures)) for i, temp in enumerate(unique_temperatures)}

# Loop over selected misorientation angles
for angle in selected_angles:
    subset = df[np.isclose(df["Misorientation Angle [°]"], angle)]

    fig, ax = plt.subplots(figsize=(8, 6))

    for temp in unique_temperatures:
        temp_subset = subset[subset["Temperature [K]"] == temp]
        temp_subset = temp_subset.sort_values(by="Stress [GPa]")

        if not temp_subset.empty:
            ax.plot(temp_subset["Stress [GPa]"],
                    temp_subset["Time Average Velocity [Ang/ps]"]*100,
                    label=f"{temp} K",
                    color=color_map[temp],
                    linestyle='--', marker='o', markersize=6)

    ax.set_xlabel("Resolved Shear Stress [GPa]", fontsize=16)
    ax.set_ylabel("Time Average Velocity [m/s]", fontsize=16)
    ax.set_title(f"Misorientation = {angle:.1f}°", fontsize=16)
    ax.legend(title="Temperature [K]", fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)

    # Save figure
    plot_filename = f"Velocity_vs_Stress_Misorientation_{angle:.1f}deg.pdf"
    plt.savefig(plot_filename, dpi=1000)
    print(f"Saved {plot_filename}")

plt.show()
