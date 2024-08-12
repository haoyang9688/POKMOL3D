import os

def process_diversity_files(base_dir, diversity_type, output_file, mean_output_file):
    all_values = []

    for i in range(32):
        folder_path = os.path.join(base_dir, str(i))

        if not os.path.exists(folder_path):
            #print(f" {folder_path} non existent. skip")
            continue

        if diversity_type == 'morgan':
            file_name = 'morgan-interal-diversity.txt'
        else:
            file_name = 'scaffold-interal-diversity.txt'
        
        file_path = os.path.join(folder_path, file_name)
        
        if os.path.exists(file_path):  
            with open(file_path, 'r') as infile:
                for line in infile:
                    if ':' in line:
                
                        value = line.split(':', 1)[1].strip()
                        all_values.append(value)

    
    if not all_values:
        print(f" ")
    else:
       
        with open(output_file, 'w') as outfile:
            for value in all_values:
                outfile.write(value + '\n')
        
        
        numeric_values = [float(value) for value in all_values if value.replace('.', '', 1).isdigit()]
        mean_value = sum(numeric_values) / len(numeric_values) if numeric_values else 0
        
        with open(mean_output_file, 'w') as mean_file:
            mean_file.write(f"Mean {diversity_type} Internal Diversity:\t{mean_value}\n")

def calculate_tdiv_merge(config):
    
    base_morgan_dir = os.path.join(config['output']['output_path'], 'General-molecular-quality-metrics', 'Morgan-Interal-diversity')
    morgan_output_file = os.path.join(base_morgan_dir, 'morgan-tdiv.txt')
    morgan_mean_output_file = os.path.join(base_morgan_dir, 'morgan-tdiv-mean.txt')
    
    base_scaffold_dir = os.path.join(config['output']['output_path'], 'General-molecular-quality-metrics', 'Scaffold-Interal-diversity')
    scaffold_output_file = os.path.join(base_scaffold_dir, 'scaffold-tdiv.txt')
    scaffold_mean_output_file = os.path.join(base_scaffold_dir, 'scaffold-tdiv-mean.txt')

    
    process_diversity_files(base_morgan_dir, 'morgan', morgan_output_file, morgan_mean_output_file)
    process_diversity_files(base_scaffold_dir, 'scaffold', scaffold_output_file, scaffold_mean_output_file)
