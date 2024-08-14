import os
import re
import pandas as pd
import numpy as np

def extract_affinity_from_txt(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if "Affinity:" in line:
                match = re.search(r'Affinity:\s*([-\d.]+)', line)
                if match:
                    return float(match.group(1))
    return None

def merge_in_situ_docking_scores(base_dir, output_file):
    docking_scores_dict = {}
    
    for i in range(32):
        subfolder_name = str(i)
        subfolder_path = os.path.join(base_dir, subfolder_name, 'qvina-text')
        
        if os.path.exists(subfolder_path):
            docking_scores = []
            for file_name in os.listdir(subfolder_path):
                if file_name.endswith('.txt'):
                    file_path = os.path.join(subfolder_path, file_name)
                    affinity = extract_affinity_from_txt(file_path)
                    if affinity is not None:
                        docking_scores.append(affinity)
            docking_scores_dict[subfolder_name] = docking_scores
        else:
            pass

    with open(output_file, 'w') as output_file:
        output_file.write('\t'.join([str(i) for i in range(32)]) + '\n')
        
        max_length = max(len(scores) for scores in docking_scores_dict.values())
        
        for i in range(max_length):
            data_row = [
                str(docking_scores_dict.get(str(j), [''])[i]) if i < len(docking_scores_dict.get(str(j), [''])) else ''
                for j in range(32)
            ]
            output_file.write('\t'.join(data_row) + '\n')

def calculate_statistics(input_file, output_dir):
    df = pd.read_csv(input_file, sep='\t', header=0)
    
    percentages = (df[df < 0].count() / df.count() * 100).fillna(0)
    percentage_output_file_path = os.path.join(output_dir, 'vina-fraction-of-PM.txt')
    with open(percentage_output_file_path, 'w') as percentage_output_file:
        percentage_output_file.write('Targets\tpositive-molecules(%)\n')
        for target, percentage in percentages.items():
            percentage_output_file.write(f'{target}\t{percentage:.3f}\n')

    df_filtered = df[df.columns[0:32]].apply(lambda x: x[x <= 0])
    df_filtered = df_filtered.dropna(axis=1, how='all')
    df_reordered = pd.DataFrame({col: df_filtered[col].dropna().reset_index(drop=True) for col in df_filtered.columns})
    means = df_reordered.mean()

    mean_output_file_path = os.path.join(output_dir, 'vina-in-situ-docking-score-PM.txt')
    with open(mean_output_file_path, 'w') as mean_output_file:
        mean_output_file.write('Targets\tMean-value\n')
        for target, mean in means.items():
            mean_output_file.write(f'{target}\t{mean:.3f}\n')
            
    def top_10_min_negative_values_mean(series):
        negative_values = series[series < 0]
        top_10_values = negative_values.nsmallest(10)
        return top_10_values.mean() if not top_10_values.empty else np.nan

    top10_means_dict = {col: top_10_min_negative_values_mean(df_reordered[col]) for col in df_reordered.columns}
    top10_output_file_path = os.path.join(output_dir, 'vina-in-situ-docking-score-top-10.txt')
    with open(top10_output_file_path, 'w') as top10_output_file:
        top10_output_file.write('Targets\tTop-10-values(mean)\n')
        for target, mean_value in top10_means_dict.items():
            top10_output_file.write(f'{target}\t{mean_value:.3f}\n')

def calculate_score(config):
    config_output_path = config['output']['output_path']
    base_dir = os.path.join(config_output_path, "Target-binding-metrics", "Vina-In-situ-Docking-results")
    output_file = os.path.join(base_dir, 'docking_scores.txt')
    out_txt_path = os.path.join(config_output_path, "Target-binding-metrics")
    merge_in_situ_docking_scores(base_dir, output_file)
    calculate_statistics(output_file, out_txt_path)
