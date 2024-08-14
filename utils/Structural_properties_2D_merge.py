import os

def calculate_structural_properties_2D_merge(config):
    config_output_path = config['output']['output_path']
    source_dir = os.path.join(config_output_path, "Structural-properties-metrics", "2D")

    target_file = os.path.join(source_dir, "2D-merge.txt")

    with open(target_file, 'w') as outfile:  
        for i in range(32):  
            folder_name = os.path.join(source_dir, str(i))  
            if os.path.exists(folder_name):  
                for file_name in os.listdir(folder_name): 
                    if file_name.endswith('.txt'):  
                        file_path = os.path.join(folder_name, file_name)  
                        with open(file_path, 'r') as infile: 
                            outfile.write(infile.read())  
            else:
                print(f"Folder {folder_name} does not exist, skipping.")

