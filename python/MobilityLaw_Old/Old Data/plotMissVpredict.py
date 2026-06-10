import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

# Directory containing the CSV files
directory = '/Users/matthewmaron/Documents/MoDELib2/python/MobilityLaw/PredictionData3'
pattern = os.path.join(directory, '*.csv')

data_dict = {}

# Define color mapping for each temperature
color_mapping = {
    100: 'royalblue',
    300: 'orange',
    500: 'forestgreen',
    700: 'firebrick',
    1000: 'purple'
}

# Read all CSV files and process them
for file_path in glob.glob(pattern):
    filename = os.path.basename(file_path)
    try:
        stress, temp = filename.replace('.csv', '').split('_')
        stress = float(stress.replace('GPa', ''))
        temp = int(temp.replace('K', ''))
    except ValueError:
        print(f"Filename format incorrect for {filename}, skipping...")
        continue
    
    df = pd.read_csv(file_path)
    
    # Ensure 'Predicted Velocity (A/ps)' contains scalar values
    if df['Predicted Velocity (A/ps)'].dtype == object:
        df['Predicted Velocity (A/ps)'] = df['Predicted Velocity (A/ps)'].apply(lambda x: x.strip('[]').split())
        df['Predicted Velocity (A/ps)'] = df['Predicted Velocity (A/ps)'].apply(lambda x: float(x[0]) if isinstance(x, list) and len(x) == 1 else float(x))
    df = df.sort_values(by=['Misorientation Angle (degrees)', 'Predicted Velocity (A/ps)'])
    
    # Initialize dictionary for the stress if not already present
    if stress not in data_dict:
        data_dict[stress] = {}

    # Store DataFrame in the dictionary
    data_dict[stress][temp] = df

# Plotting
for stress, temp_data in data_dict.items():
    plt.figure(figsize=(8, 6))
    
    for temp in sorted(temp_data.keys()):  # Sort by temperature to ensure order in legend
        df = temp_data[temp]
        color = color_mapping.get(temp, 'black')  # Default to black if temp is not in color_mapping
        plt.scatter(df['Misorientation Angle (degrees)'], df['Predicted Velocity (A/ps)'], label=f'{temp}K', color=color)
        plt.plot(df['Misorientation Angle (degrees)'], df['Predicted Velocity (A/ps)'], linestyle='--', color=color)

    plt.title(f'Applied Stress {stress} GPa')
    plt.xlabel('Misorientation Angle (degrees)')
    plt.ylabel('Predicted Velocity (A/ps)')
    plt.ylim(0,1.75)
    plt.xlim(0,180)
    plt.legend(title='Temperature', loc='best')
    plt.grid(True)
    plt.savefig(f'PredictionDataFigs3/Misorientation_vs_Velocity_Stress_{stress}GPa.png', format='png')
    plt.show()
