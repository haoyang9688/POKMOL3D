import os
import gzip
import shutil
import numpy as np
import pandas as pd

def calculate_in_situ_docking_score(config):
    config_output_path = config['output']['output_path']
    res_omaegz_root = os.path.join(config_output_path, "Target-binding-metrics", "Glide-In-situ-Docking-results")
    out_txt_path = os.path.join(config_output_path, "Target-binding-metrics")
    
    #Decompression and extraction of docking score
    for number_folder in os.listdir(res_omaegz_root):
        number_folder_path = os.path.join(res_omaegz_root, number_folder)
        docking_txt_folder = os.path.join(res_omaegz_root, number_folder, "In-situ-Docking-score")
        
        if os.path.exists(number_folder_path) and os.path.isdir(number_folder_path):
            
            docking_txt_path = os.path.join(docking_txt_folder, "in-situ-docking-score.txt")
            if not os.path.exists(docking_txt_folder):
                os.makedirs(docking_txt_folder)
            
            #Traverse all. sdfgz files under numb_folder_path
            for sdfgz_file_name in os.listdir(number_folder_path):
                if sdfgz_file_name.endswith(".sdfgz"):
                    sdfgz_file_path = os.path.join(number_folder_path, sdfgz_file_name)
                    
                    #Extract sdfgz file
                    with gzip.open(sdfgz_file_path, 'rb') as f_in:
                        sdf_file_name = os.path.splitext(sdfgz_file_name)[0] + ".sdf"
                        sdf_file_path = os.path.join(number_folder_path, sdf_file_name)
                        with open(sdf_file_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    
                    #Read the decompressed sdf file and extract<r_i_docking_Score>
                    with open(sdf_file_path, 'r') as sdf_file:
                        lines = sdf_file.readlines()
                        docking_score = None
                        for i in range(len(lines)):
                            if lines[i].strip() == "> <r_i_docking_score>":
                                docking_score = lines[i+1].strip()
                                break
                    
                    #Write the docking-score.txt file
                    if docking_score:
                        sdf_prefix = os.path.splitext(sdfgz_file_name)[0]  # 取文件名前缀
                        with open(docking_txt_path, 'a') as docking_txt_file:
                            docking_txt_file.write(f"{sdf_prefix}:{docking_score}\n")

    docking_scores_dict = {}
    for i in range(32):
        subfolder_name = str(i)
        subfolder_path = os.path.join(res_omaegz_root, subfolder_name, 'In-situ-Docking-score')
        if os.path.exists(subfolder_path):
            docking_scores_file_path = os.path.join(subfolder_path, 'in-situ-docking-score.txt')
            if os.path.exists(docking_scores_file_path):
                with open(docking_scores_file_path, 'r') as file:
                    docking_scores = [float(line.strip().split(':')[-1]) for line in file]
                    docking_scores_dict[subfolder_name] = docking_scores
    output_file_path = os.path.join(res_omaegz_root, 'in-situ-docking-score-merge.txt')
    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join([str(i) for i in range(32)]) + '\n')

        max_length = max(len(scores) for scores in docking_scores_dict.values())
        
        #Write data
        for i in range(max_length):
            #Extract the i-th score for each folder, and replace it with an empty string if the number of scores in the folder is insufficient
            data_row = [str(docking_scores_dict.get(str(j), [''])[i]) if i < len(docking_scores_dict.get(str(j), [''])) else '' for j in range(32)]
            output_file.write('\t'.join(data_row) + '\n')

    #Read the docking-scoremerging. txt file and count the proportion of values less than 0 in each column
    df = pd.read_csv(output_file_path, sep='\t')
    percentages = (df[df < 0].count() / df.count() * 100).fillna(0)

    #Output the proportion of values less than 0 in each column as a txt file
    percentage_output_file_path = os.path.join(out_txt_path, 'glide-fraction-of-PM.txt')
    with open(percentage_output_file_path, 'w') as percentage_output_file:
        percentage_output_file.write('Targets\tpositive-molecules(%)\n')
        for target, percentage in percentages.items():
            percentage_output_file.write(f'{target}\t{percentage:.3f}\n')

    #Extract non 10000 values and exclude positive values
    df_filtered = df[df.columns[0:32]].apply(lambda x: x[(x != 10000) & (x <= 0)])

    #Remove null values
    df_filtered = df_filtered.dropna(axis=1, how='all')

    #Rearrange the values of each column
    df_reordered = pd.DataFrame({
        col: df_filtered[col].dropna().reset_index(drop=True)
        for col in df_filtered.columns
    })

    #Calculate the average value of each column
    means = df_reordered.mean()

    #Output the average result as a txt file
    mean_output_file_path = os.path.join(out_txt_path, 'glide-in-situ-docking-score-PM.txt')
    with open(mean_output_file_path, 'w') as mean_output_file:
        mean_output_file.write('Targets\tMean-value\n')
        for target, mean in means.items():
            mean_output_file.write(f'{target}\t{mean:.3f}\n')

    #Define a function to extract the top 10 smallest negative values of each column and calculate the average value
    def top_10_min_negative_values_mean(series):
        negative_values = series[series < 0]
        top_10_values = negative_values.nsmallest(10)
        return top_10_values.mean() if not top_10_values.empty else np.nan

    #Extract the average of the top 10 smallest negative values in each column
    top10_means_dict = {col: top_10_min_negative_values_mean(df_reordered[col]) for col in df_reordered.columns}

    #Output the average of the first 10 smallest negative values as a txt file
    top10_output_file_path = os.path.join(out_txt_path, 'glide-in-situ-docking-score-top-10.txt')
    with open(top10_output_file_path, 'w') as top10_output_file:
        top10_output_file.write('Targets\tTop-10-values(mean)\n')
        for target, mean_value in top10_means_dict.items():
            top10_output_file.write(f'{target}\t{mean_value:.3f}\n')

