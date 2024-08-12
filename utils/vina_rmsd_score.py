import os
import pandas as pd

def calculate_average_from_file(file_path):
    """
    Calculate the average value of the numerical data in the given file, ignoring lines with a value of -1。
    
    :param file_path: File Path
    :return: mean value
    """
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
    """
    Extract the top 10 non negative minimum values from each column and calculate their average。
    
    :param series: Data Column
    :return: The average of the top 10 non negative minimum values
    """
    #Filter out negative values
    positive_values = series[series >= 0]
    top10_values = positive_values.nsmallest(10)
    if not top10_values.empty:
        return top10_values.mean()
    else:
        return float('nan')

def glide_rmsd_score_summary(config):

    config_output_path = config['output']['output_path']
    rmsd_method = config['settings']['RMSD_settings']['rmsd_method']
    if rmsd_method == "openeye":
        root_folder_path = os.path.join(config_output_path, "Target-binding-metrics", "Vina-Openeye-RMSD-results")
    if rmsd_method == "rdkit":
        root_folder_path = os.path.join(config_output_path, "Target-binding-metrics", "Vina-RDKit-RMSD-results")
    output_txt_path = os.path.join(config_output_path, "Target-binding-metrics")
    os.makedirs(root_folder_path, exist_ok=True) 
  
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

    #print(f"extracted docking Score saved in {output_file_path}")

    
    average_output_file_path = os.path.join(output_txt_path, 'rmsd-score.txt')
    with open(average_output_file_path, 'w') as average_output_file:
        average_output_file.write("Targets\tMean-value\n")
        for col in range(32):
            col_values = docking_scores_dict.get(str(col), [])
            
            col_values = [float(v) for v in col_values if v not in ['', '-1']]
            avg_value = sum(col_values) / len(col_values) if col_values else float('nan')
            average_output_file.write(f"{col}\t{avg_value:.3f}\n")

    #print(f"mean value saved in {average_output_file_path}")

   
    df = pd.read_csv(output_file_path, sep='\t', header=0)  

    
    top10_means_dict = {col: top_10_min_positive_values_mean(df[col]) for col in df.columns}

    
    top10_means_df = pd.DataFrame(list(top10_means_dict.items()), columns=['Targets', 'Top-10-values(mean)'])
    top10_output_file_path = os.path.join(output_txt_path, 'rmsd-score-top10.txt')
    top10_means_df.to_csv(top10_output_file_path, sep='\t', header=True, index=None)

    #print(f"The average of the top 10 non negative minimum values saved to {top10_output_file_path}")
