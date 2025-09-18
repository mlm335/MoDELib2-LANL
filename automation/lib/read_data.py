import subprocess
import os
import pathlib
import sys
import numpy as np
import string, os, math, sys

# Copy folder contents to a new folder
def copy_folder_contents(src_folder, dest_folder):
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
    
    command = ['cp', '-r', os.path.join(src_folder, '.'), dest_folder]
    subprocess.run(command, check=True)
    print(f"Copied contents of {src_folder} to {dest_folder}")


# Generate a new folder in a directory
def generate_folder(base_dir, new_folder_name):
    new_folder_path = os.path.join(base_dir, new_folder_name)
    
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        
    return new_folder_path
    
    
#Generate file path from directory and file name
def generate_file_paths(directory, names):
    file_paths = []
    for name in names:
        file_path = os.path.join(directory, name)
        file_paths.append(file_path)
    return file_paths


# Get variable from a txt file
def get_scalar(file_path, variable_name):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.split(';', 1)[0]  # Ignore everything after the first semicolon
                parts = line.strip().split('=')
                if len(parts) == 2:
                    name, value = parts
                    if name.strip() == variable_name:
                        print(f"Found variable {variable_name} in the line: {line.strip()}")
                        return float(value.strip())
        return None  # Variable not found in the file
    except FileNotFoundError:
        return None  # File not found

# Get vector from a txt file
def get_vector(file_path, variable_name):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.split(';', 1)[0]  # Ignore everything after the first semicolon
                parts = line.strip().split('=')
                if len(parts) == 2:
                    name, value = parts
                    if name.strip() == variable_name:
                        print(f"Found variable {variable_name} in the line: {line.strip()}")
                        try:
                            vector_str = value.strip()
                            # Add commas between values
                            vector_str_with_commas = ','.join(vector_str.split())
                            vector = np.array([float(item) for item in vector_str_with_commas.split(',')])
                            return vector
                        except ValueError as e:
                            print(f"Error converting to a NumPy array: {e}")
                            return None
        return None  # Variable not found in the file
    except FileNotFoundError:
        return None  # File not found


#def get_vector(file_path, variable_name):
#    try:
#        with open(file_path, 'r') as file:
#            for line in file:
#                line = line.split(';', 1)[0]  # Ignore everything after the first semicolon
#                parts = line.strip().split('=')
#                if len(parts) == 2:
#                    name, value = parts
#                    if name.strip() == variable_name:
#                        print(f"Found variable {variable_name} in the line: {line.strip()}")
#                        try:
#                            vector_str = value.strip()
#                            vector = [float(item) for item in vector_str.split(',')]
#                            return vector
#                        except ValueError:
#                            return None  # Unable to convert to a vector
#        return None  # Variable not found in the file
#    except FileNotFoundError:
#        return None  # File not found

# Read F File
def readFfile(folder):
    F=np.loadtxt(folder +'/F/F_0.txt');
    with open(folder +'/F/F_labels.txt') as f:
        lines = f.readlines()
        for idx in range(len(lines)):
            lines[idx] = lines[idx].rstrip()
    return F,lines

# Get array from F file
def getFarray(F,Flabels,label):
    k=0;
    for line in Flabels:
        if line==label:
            return F[:,k]
        k=k+1
    return np.zeros(shape=(0,0))

def load_variables_from_file(filename):
    variables = {}
    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line and '=' in line:
                     semicolon_index = line.find(';')
                     if semicolon_index != -1:
                    # Cut the line at the first semicolon
                        line = line[:semicolon_index] + '\n'
                     name, value = line.split('=', 1)
                     #value = value.replace('e', 'e')
                     #value = value.replace('-', '-')
                     #value = value.replace('+', '+0')
                     variables[name.strip()] = value.strip()
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return variables

def copy_file_with_variable_changes(source_file, destination_file, changes):
    try:
        with open(source_file, 'r') as src, open(destination_file, 'w') as dest:
            for line in src:
                for var_name, new_value in changes.items():
                    line = line.replace(f'{var_name} =', f'{var_name} = {new_value}')
                dest.write(line)

        print("File copied with variable changes.")
    except FileNotFoundError:
        print(f"Source file '{source_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
