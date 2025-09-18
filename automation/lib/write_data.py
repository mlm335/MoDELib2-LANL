def modify_file(file_path, modifications):
    # Read the content of the polycrystal file
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Iterate through modifications and apply changes
    for key, value in modifications.items():
        key_index = None

        # Search for the key in the polycrystal file
        for idx, line in enumerate(lines):
            if key in line:
                key_index = idx
                break
        
        # If the key is found, modify the lines
        if key_index is not None:
            if isinstance(value, list):  # For matrices like F
                # Flatten the matrix and replace the next lines
                for i in range(len(value)):
                    formatted_line = ' '.join([' '.join(sublist) for sublist in value[i]]) + '\n'
                    lines[key_index + 1 + i] = formatted_line
            else:  # For single line entries like meshFile
                lines[key_index] = f'{key}: {value}\n'

    # Write the modified content back to the polycrystal file
    with open(file_path, 'w') as f:
        f.writelines(lines)
        

def format_value(value):
    """
    Formats the value based on its inferred type (scalar, vector, or matrix).

    Args:
    - value (any): The value to format.

    Returns:
    - str: The formatted value as a string.
    """
    if isinstance(value, list):
        if all(isinstance(i, list) for i in value):  # Matrix
            return '\n'.join([' '.join(map(str, row)) for row in value]) + ';'
        else:  # Vector
            return ' '.join(map(str, value)) + ';'
    else:  # Scalar
        return f"{value};"
