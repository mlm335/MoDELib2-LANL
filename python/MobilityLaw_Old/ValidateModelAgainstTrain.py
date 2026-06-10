import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
import pandas as pd
from matplotlib import rc
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['font.family'] = ['Times New Roman']
plt.rcParams['text.usetex'] = True
rc('text', usetex=True)
rc('font', family='serif')

# -------------------------------
# Load FULL trained pipeline
# -------------------------------
loaded_model = joblib.load("GBR_model.pkl")

# -------------------------------
# Set up Simulation
# -------------------------------
df = pd.read_excel('aug20_use_for_velocity.xlsx')

target_e1 = [0.1, 0.25, 0.4, 0.5, 0.75]
target_T = [100, 300, 500, 700, 1000]
target_miss= [0,11,31,48,67,90,110,120,152,174]

# Nested dict: {miss: {temp: [(e1, velocity), ...]}}
results_by_miss = {miss: {T: [] for T in target_T} for miss in target_miss}

for T in target_T:
    df_T = df[df['temp'] == T]  # Filter rows for temperature

    for miss in target_miss:
        df_T_miss=df_T[df_T['miss'] == miss] # Filter rows for miss angle

        for stress in target_e1:
            if not df_T_miss.empty:
                idx = (df_T_miss['e1'] - stress).abs().idxmin()
                row = df_T_miss.loc[idx]

                print(f"Found Data For T={T}, miss={miss}, and e1={stress} in row {idx+2}")
                
                features = pd.DataFrame([{
                    'e1': row['e1'],
                    'e2': row['e2'],
                    'e3': row['e3'],
                    'anev1': row['anev1'],
                    'anev2': row['anev2'],
                    'anev3': row['anev3'],
                    'temp': row['temp'],
                    'miss': row['miss']
                }])

                predicted_velocity = loaded_model.predict(features)
                measured_velocity = row['velo']  

                results_by_miss[miss][T].append((row['e1'], predicted_velocity*100, measured_velocity*100))
            else:
                print(f"No match for temp={T}, miss={miss}")


# -------------------------------
# Plot Results
# -------------------------------
cmap = plt.get_cmap("viridis")
color_map = {temp: cmap(i / len(target_T)) for i, temp in enumerate(target_T)}
temperature_colors = {
    100: '#1f77b4',  # Blue
    300: '#ff7f0e',  # Orange
    500: '#2ca02c',  # Green
    700: '#d62728',  # Red
    1000: '#9467bd'  # Purple
}

for miss in target_miss:
    plt.figure(figsize=(8, 6))
    for T in target_T:
        data = results_by_miss[miss][T]
        if data:
            stresses, predicted_v, measured_v = zip(*data)
            plt.plot(stresses, predicted_v, 'x--',color=temperature_colors[T], label=f'Predicted {T}K')
            plt.plot(stresses, measured_v, 'o-',color=temperature_colors[T], label=f'Measured {T}K')
    
    plt.xlabel(r'$e_1$ Value [GPa]', fontsize=20)
    plt.ylabel(r'Velocity [m/s]', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tick_params(axis='both', which='minor', labelsize=16)
    plt.xlim(0,1)
    plt.title(rf'Missorientation = {miss}$^\circ$', fontsize=24)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'New_Validation_{miss}degrees_ValidateAgainstTrain.pdf', dpi=1000)
plt.show()
