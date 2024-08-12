# utils/SA_merge.py
import os

def calculate_SA_merge(config):
    
    source_dir = os.path.join(config['output']['output_path'], 'General-molecular-quality-metrics', 'SA')

    target_dir = os.path.join(config['output']['output_path'], 'General-molecular-quality-metrics', 'SA')

    column_values = {'SA Score:': []}

   
    for i in range(32): 
        folder_path = os.path.join(source_dir, str(i))
        if not os.path.exists(folder_path):
            print(f"{folder_path} non-existent. skip")
            continue

        for file_name in os.listdir(folder_path):
            if file_name.endswith('.txt'): 
                file_path = os.path.join(folder_path, file_name)
                with open(file_path, 'r') as infile:
                       
                    for line in infile:
                           
                        if 'SA Score:' in line:
                                
                            score = line.split('SA Score:')[1].strip()
                            column_values['SA Score:'].append(score)

    
    target_file_path = os.path.join(target_dir, "all-SA-Score.txt")
    with open(target_file_path, 'w') as outfile:
        for value in column_values['SA Score:']:
            outfile.write(value + '\n')
    
    scores = [float(score) for score in column_values['SA Score:']]
    mean_value = sum(scores) / len(scores) if scores else 0
    mean_file_path = os.path.join(target_dir, 'SA-mean.txt')
    with open(mean_file_path, 'w') as mean_file:
        mean_file.write(f"Mean SA Score: {mean_value}\n")
