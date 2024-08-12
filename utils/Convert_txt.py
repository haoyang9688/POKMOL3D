import os

def convert_txt(folder_path):
    #Check if the folder exists and is a directory
    if os.path.isdir(folder_path):
        for filename in os.listdir(folder_path):
            if filename.endswith('.smi'):
                input_file_path = os.path.join(folder_path, filename)
                output_file_path = os.path.join(folder_path, os.path.splitext(filename)[0] + '.txt')

                #Open the original file for reading and write a new file
                with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
                    #Skip the first line header
                    next(input_file)

                    #Read each line, convert the format and write it to a new file
                    for line in input_file:
                        line = line.strip()
                        if line:
                            split_line = line.split()
                            if len(split_line) == 2:
                                name, smiles = split_line
                                output_file.write(f"{smiles}\t{name}\n")
    else:
        print(f"Directory {folder_path} not found.")

def calculate_convert_txt(config):
    config_output_path = config['output']['output_path']
    input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics')  
    uniqueness_path = os.path.join(input_path, 'Uniqueness')
    for folder_number in range(32):  
        folder_path = os.path.join(uniqueness_path, str(folder_number))
        if os.path.isdir(folder_path):
            convert_txt(folder_path)
        else:
            print(f"Subdirectory {folder_path} not found. Skipping...")
