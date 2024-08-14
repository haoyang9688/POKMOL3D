import os

def calculate_average_scores(base_dir, average_output_file):
    results = []

    for i in range(32):
        folder_path = os.path.join(base_dir, str(i))
        file_path = os.path.join(folder_path, 'highest_score.txt')

        if os.path.exists(file_path):
            scores = []
            with open(file_path, 'r') as file:
                for line in file:
                    parts = line.strip().split()
                    if parts and parts[-1].replace('.', '', 1).replace('-', '', 1).isdigit():  #Check if the last column is a number
                        score = float(parts[-1])
                        scores.append(score)
            if scores:
                mean_value = sum(scores) / len(scores)
                results.append((i, mean_value))
        else:
            pass

    with open(average_output_file, 'w') as outfile:
        outfile.write("Targets\tmean-value\n")
        for target, mean_value in results:
            outfile.write(f"{target}\t{mean_value:.3f}\n")

    #print(f"The average value has been calculated and written in {average_output_file}")

def calculate_top_10_average_scores(base_dir, top_10_output_file):
    results = []

    for i in range(32):
        folder_path = os.path.join(base_dir, str(i))
        file_path = os.path.join(folder_path, 'highest_score.txt')

        if os.path.exists(file_path):
            scores = []
            with open(file_path, 'r') as file:
                for line in file:
                    parts = line.strip().split()
                    if parts and parts[-1].replace('.', '', 1).replace('-', '', 1).isdigit():  #Check if the last column is a number
                        score = float(parts[-1])
                        if score < 0:  #Ensure that only negative values are processed
                            scores.append(score)
            if scores:
                #Take the first 10 smallest negative values. If there are less than 10, use the actual quantity
                top_10_scores = sorted(scores)[:10]
                mean_top_10_value = sum(top_10_scores) / len(top_10_scores)
                results.append((i, mean_top_10_value))
        else:
            pass
    with open(top_10_output_file, 'w') as outfile:
        outfile.write("Targets\tTop-10-values(mean)\n")
        for target, mean_value in results:
            outfile.write(f"{target}\t{mean_value:.3f}\n")

    #print(f"The average value of Top-10 has been calculated and written in {top_10_output_file}")

def calculate_docking_score(config):
    base_dir = os.path.join(config['output']['output_path'], 'Target-binding-metrics', 'Vina-Docking-results')
    average_output_file = os.path.join(config['output']['output_path'], 'Target-binding-metrics', 'vina-docking-score.txt')
    top_10_output_file = os.path.join(config['output']['output_path'], 'Target-binding-metrics', 'vina-docking-score-top10.txt')
    os.makedirs(config['output']['output_path'], exist_ok=True)
    
    calculate_average_scores(base_dir, average_output_file)
    calculate_top_10_average_scores(base_dir, top_10_output_file)
