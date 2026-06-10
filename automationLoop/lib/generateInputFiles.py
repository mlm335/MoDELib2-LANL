def writeInputFiles(file_name, key, datatochange):
    print(f"Changing {key} to {datatochange} in {file_name}")
    
    with open(file_name, 'r') as f:
        lines = f.readlines()
        
        for i, line in enumerate(lines):
            if line.strip().startswith(f"{key}="):  # If the line starts with the key
                
                if isinstance(datatochange, list) and all(isinstance(row, list) for row in datatochange):
                    comment = ''
                    if ';' in lines[i + len(datatochange) - 1] and '#' in lines[i + len(datatochange) - 1]:
                        comment = lines[i + len(datatochange) - 1].split('#', 1)[1].strip()
                 
                    lines[i] = f"{key}={' '.join(datatochange[0])}\n"
                    for j in range(1, len(datatochange) - 1):
                        lines[i + j] = ' '.join(datatochange[j]) + '\n'
                    lines[i + len(datatochange) - 1] = ' '.join(datatochange[-1]) + f"; {'# ' + comment if comment else ''}\n"
                
                elif isinstance(datatochange, list) and all(isinstance(item, str) for item in datatochange):
                    comment = ''
                    if ';' in line and '#' in line:
                        comment = line.split('#', 1)[1].strip()
                    lines[i] = f"{key}={' '.join(datatochange)}; {'# ' + comment if comment else ''}\n"

                elif isinstance(datatochange, str):
                    comment = ''
                    if ';' in line and '#' in line:
                        comment = line.split('#', 1)[1].strip()  
                    lines[i] = f"{key}={datatochange}; {'# ' + comment if comment else ''}\n"

                else:
                    print("Data type could not be determined, change input file")
                    
                break

    with open(file_name, 'w') as f:
        f.writelines(lines)
