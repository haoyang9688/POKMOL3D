import os

def calculate_uniqueness(config):
 
    config_output_path = config['output']['output_path']
    input_folder = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Validity')
    output_base_folder = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness')

    
    os.makedirs(output_base_folder, exist_ok=True)

    for dir_number in range(32):
        dir_path = os.path.join(input_folder, str(dir_number))
        output_folder = os.path.join(output_base_folder, str(dir_number))
        
    
        if not os.path.isdir(dir_path):
            #print(f"Directory {dir_path} not found. Skipping...")
            continue
        
       
        for root, _, files in os.walk(dir_path):
            for file in files:
                if file.endswith(".smi"):
                    input_file = os.path.join(root, file)
                    output_file = os.path.join(output_folder, "unique.smi")

                    
                    os.makedirs(output_folder, exist_ok=True)
                    
                    with open(input_file, 'r') as f:
                        lines = f.readlines()
                    
                    molmap = set()
                    new_lines = []
                    valcount = 0
                    uniqcount = 0
                    
                    for line in lines:
                        valcount += 1
                        smiles_str = line.split(" ")[0]
                        if smiles_str not in molmap:
                            new_lines.append(line)
                            molmap.add(smiles_str)
                            uniqcount += 1

                    with open(output_file, 'w') as uniqoutputfile:
                        for line in new_lines:
                            uniqoutputfile.write(line + '\n')

                    #print(f"Processed {input_file}. Output saved to {output_file}")