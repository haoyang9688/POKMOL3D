import os

def convert_smi(config):
    """
    Traverse subfolders 0-31 along the path,
Convert the. txt files in each subfolders to. smi files.
    """
    
    output_path = config['output']['output_path']
    input_folder = os.path.join(output_path, 'General-molecular-quality-metrics','Validity')
    #Traverse folders 0-31
    for dir_number in range(32):
        dir_name = str(dir_number)
        dir_path = os.path.join(input_folder, dir_name)

        # Check if the folder exists and is a directory
        if not os.path.isdir(dir_path):
            #print(f"Directory {dir_path} not found. Skipping...")
            continue

        for filename in os.listdir(dir_path):
            if filename.endswith('.txt'):
                input_file_path = os.path.join(dir_path, filename)
                output_file_path = os.path.join(dir_path, os.path.splitext(filename)[0] + '.smi')

                #Open the original file for reading and write a new file
                with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
                    #Write the header of the new file
                    output_file.write("SMILES Name\n")

                    #Read each line, convert the format and write it to a new file
                    for line in input_file:
                        line = line.strip()
                        if line:
                            split_line = line.split('\t')
                            if len(split_line) == 2:
                                name, smiles = split_line
                                output_file.write(f"{smiles} {name}\n")
    else:
        print(f"Directory {dir_path} not found.")
