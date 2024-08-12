import os
import gzip
import shutil
import numpy as np
import pandas as pd

def calculate_docking_score(config):
    # Define root path
    config_output_path = config['output']['output_path']
    res_omaegz_root = os.path.join(config_output_path, "Target-binding-metrics", "Glide-Docking-results")
    score_txt_path = os.path.join(config_output_path, "Target-binding-metrics")
    
    # Decompression and extraction of docking score
    for number_folder in os.listdir(res_omaegz_root):
        number_folder_path = os.path.join(res_omaegz_root, number_folder)
        docking_txt_folder = os.path.join(res_omaegz_root, number_folder, "Docking-score")
        
        if os.path.exists(number_folder_path) and os.path.isdir(number_folder_path):
            
            docking_txt_path = os.path.join(docking_txt_folder, "docking-score.txt")
            if not os.path.exists(docking_txt_folder):
                os.makedirs(docking_txt_folder)
            
            for sdfgz_file_name in os.listdir(number_folder_path):
                if sdfgz_file_name.endswith(".sdfgz"):
                    sdfgz_file_path = os.path.join(number_folder_path, sdfgz_file_name)
                    
                    # Extract sdfgz file
                    with gzip.open(sdfgz_file_path, 'rb') as f_in:
                        sdf_file_name = os.path.splitext(sdfgz_file_name)[0] + ".sdf"
                        sdf_file_path = os.path.join(number_folder_path, sdf_file_name)
                        with open(sdf_file_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    
                    # Read the decompressed sdf file and extract <r_i_docking_Score>
                    with open(sdf_file_path, 'r') as sdf_file:
                        lines = sdf_file.readlines()
                        docking_score = None
                        for i in range(len(lines)):
                            if lines[i].strip() == "> <r_i_docking_score>":
                                docking_score = lines[i+1].strip()
                                break
                    
                    if docking_score:
                        sdf_prefix = os.path.splitext(sdfgz_file_name)[0]  # Retrieve file name prefix
                        with open(docking_txt_path, 'a') as docking_txt_file:
                            docking_txt_file.write(f"{sdf_prefix}:{docking_score}\n")

    docking_scores_dict = {}

    for i in range(32):
        subfolder_name = str(i)
        subfolder_path = os.path.join(res_omaegz_root, subfolder_name, 'Docking-score')

        if os.path.exists(subfolder_path):
            docking_scores_file_path = os.path.join(subfolder_path, 'docking-score.txt')

            if os.path.exists(docking_scores_file_path):
                with open(docking_scores_file_path, 'r') as file:
                    docking_scores = [float(line.strip().split(':')[-1]) for line in file]

                    docking_scores_dict[subfolder_name] = docking_scores
    
    output_file_path = os.path.join(res_omaegz_root, 'glide-docking-score-merge.txt')
    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join([str(i) for i in range(32)]) + '\n')

        max_length = max(len(scores) for scores in docking_scores_dict.values())
        
        for i in range(max_length):
            data_row = [str(docking_scores_dict.get(str(j), [''])[i]) if i < len(docking_scores_dict.get(str(j), [''])) else '' for j in range(32)]
            output_file.write('\t'.join(data_row) + '\n')

    # Calculate the average value of each column
    mean_values = {}
    for i in range(32):
        scores = docking_scores_dict.get(str(i), [])
        if scores:
            mean_values[i] = np.mean(scores)

    # Output the average result as a txt file
    mean_output_file_path = os.path.join(score_txt_path, 'glide-docking-score.txt')
    with open(mean_output_file_path, 'w') as mean_output_file:
        mean_output_file.write('Targets\tMean-value\n')
        for target, mean_value in mean_values.items():
            mean_output_file.write(f'{target}\t{mean_value:.3f}\n')

    df = pd.read_csv(output_file_path, sep='\t', dtype=str)
    top10_means = df.apply(lambda x: pd.to_numeric(x, errors='coerce').nsmallest(10).mean(), axis=0)

    top10_output_file_path = os.path.join(score_txt_path, 'glide-docking-score-top-10.txt')
    with open(top10_output_file_path, 'w') as top10_output_file:
        top10_output_file.write('Targets\tTop-10-values(mean)\n')
        for target, mean_value in top10_means.items():
            top10_output_file.write(f'{target}\t{mean_value:.3f}\n')

    # Print the saved file paths
    #print(f"Docking scores have been saved to {output_file_path}")
    #print(f"Mean values have been saved to {mean_output_file_path}")
    #print(f"Top 10 mean values have been saved to {top10_output_file_path}")

