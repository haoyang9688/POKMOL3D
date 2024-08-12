import os
import pandas as pd

def calculate_average_from_file(file_path):
    scores = []
    try:
        with open(file_path, 'r') as file:
            next(file)  
            for line in file:
                try:
                    score = float(line.split(':')[1].strip())
                    if score != -1:
                        scores.append(score)
                except (IndexError, ValueError) as e:
                    print(f"Error processing line in {file_path}: {line.strip()} - {e}")
        
        if scores:
            avg = sum(scores) / len(scores)
        else:
            avg = 0.0
        return avg
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None

def top_10_min_positive_values_mean(series):
    positive_values = pd.to_numeric(series, errors='coerce')
    positive_values = positive_values[positive_values >= 0]
    top10_values = positive_values.nsmallest(10)
    if not top10_values.empty:
        return top10_values.mean()
    else:
        return float('nan')

def process_rmsd_scores(root_folder_path, output_txt_path, avg_file_name, top10_file_name):
    docking_scores_dict = {}

    for i in range(32):
        subfolder_name = str(i)
        subfolder_path = os.path.join(root_folder_path, subfolder_name, 'rmsd-score')

        if os.path.exists(subfolder_path):
            docking_scores_file_path = os.path.join(subfolder_path, 'rmsd-summary.txt')

            if os.path.exists(docking_scores_file_path):
                with open(docking_scores_file_path, 'r') as file:
                    docking_scores = [line.strip().split()[-1] for line in file]
                    docking_scores_dict[subfolder_name] = docking_scores

    output_file_path = os.path.join(root_folder_path, 'rmsd-score-merge.txt')
    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join([str(i) for i in range(32)]) + '\n')
        max_length = max(len(scores) for scores in docking_scores_dict.values())
        for i in range(max_length):
            data_row = []
            for j in range(32):
                scores = docking_scores_dict.get(str(j), [])
                if i < len(scores):
                    data_row.append(scores[i])
                else:
                    data_row.append('')
            output_file.write('\t'.join(data_row) + '\n')

    df = pd.read_csv(output_file_path, sep='\t', header=0)

    average_output_file_path = os.path.join(output_txt_path, avg_file_name)
    with open(average_output_file_path, 'w') as average_output_file:
        average_output_file.write("Targets\tMean-value\n")
        for col in df.columns:
            col_values = df[col].dropna()
            col_values = col_values[col_values != -1]
            avg_value = col_values.mean() if not col_values.empty else float('nan')
            average_output_file.write(f"{col}\t{avg_value:.3f}\n")

    top10_means_dict = {col: top_10_min_positive_values_mean(df[col]) for col in df.columns}
    top10_means_df = pd.DataFrame(list(top10_means_dict.items()), columns=['Targets', 'Top-10-values(mean)'])
    top10_output_file_path = os.path.join(output_txt_path, top10_file_name)
    top10_means_df.to_csv(top10_output_file_path, sep='\t', header=True, index=None)

def rmsd_score_summary(config):
    config_output_path = config['output']['output_path']
    rmsd_method = config['settings']['RMSD_settings']['rmsd_method']
    output_txt_path = os.path.join(config_output_path, "Target-binding-metrics")
    os.makedirs(output_txt_path, exist_ok=True)
    
    if rmsd_method == "openeye":
        if os.path.exists(os.path.join(config_output_path, "Target-binding-metrics", "Glide-Openeye-RMSD-results")):
            root_folder_path = os.path.join(config_output_path, "Target-binding-metrics", "Glide-Openeye-RMSD-results")
            process_rmsd_scores(root_folder_path, output_txt_path, 'glide-openeye-rmsd-score.txt', 'glide-openeye-rmsd-score-top10.txt')
        if os.path.exists(os.path.join(config_output_path, "Target-binding-metrics", "Vina-Openeye-RMSD-results")):
            root_folder_path = os.path.join(config_output_path, "Target-binding-metrics", "Vina-Openeye-RMSD-results")
            process_rmsd_scores(root_folder_path, output_txt_path, 'vina-openeye-rmsd-score.txt', 'vina-openeye-rmsd-score-top10.txt')
    elif rmsd_method == "rdkit":
        if os.path.exists(os.path.join(config_output_path, "Target-binding-metrics", "Glide-RDKit-RMSD-results")):
            root_folder_path = os.path.join(config_output_path, "Target-binding-metrics", "Glide-RDKit-RMSD-results")
            process_rmsd_scores(root_folder_path, output_txt_path, 'glide-rdkit-rmsd-score.txt', 'glide-rdkit-rmsd-score-top10.txt')
        if os.path.exists(os.path.join(config_output_path, "Target-binding-metrics", "Vina-RDKit-RMSD-results")):
            root_folder_path = os.path.join(config_output_path, "Target-binding-metrics", "Vina-RDKit-RMSD-results")
            process_rmsd_scores(root_folder_path, output_txt_path, 'vina-rdkit-rmsd-score.txt', 'vina-rdkit-rmsd-score-top10.txt')
