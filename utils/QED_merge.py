import os

def calculate_QED_merge(config):
    source_dir = os.path.join(config['output']['output_path'], 'General-molecular-quality-metrics', 'QED')
    target_dir = os.path.join(config['output']['output_path'], 'General-molecular-quality-metrics', 'QED')

    columns_to_extract = ['QED']

    column_values = {column: [] for column in columns_to_extract}

    for i in range(32):
        folder_path = os.path.join(source_dir, str(i))
        if not os.path.exists(folder_path):
            #print(f" {folder_path} non-existent. skip。")
            continue

        for file_name in os.listdir(folder_path):
            if file_name.endswith('.txt'): 
                file_path = os.path.join(folder_path, file_name)
                with open(file_path, 'r') as infile:
        
                    header = infile.readline().strip().split('\t')
                    column_indices = {column: header.index(column) for column in columns_to_extract}
                       
                    for line in infile:
                        cols = line.strip().split('\t')
                        for column in columns_to_extract:
                            column_values[column].append(float(cols[column_indices[column]]))

    for column in columns_to_extract:
        target_file_path = os.path.join(target_dir, f"all-{column}.txt")
        with open(target_file_path, 'w') as outfile:
            for value in column_values[column]:
                outfile.write(f"{value}\n")
        
        mean_value = sum(column_values[column]) / len(column_values[column])
        mean_file_path = os.path.join(target_dir, 'QED-mean.txt')
        with open(mean_file_path, 'w') as mean_file:
            mean_file.write(f"Mean\tQED:\t{mean_value}\n")
