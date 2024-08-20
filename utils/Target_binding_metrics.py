import os
import pandas as pd

def calculate_metrics(config):
    """
    Calculate metrics for target binding.
    """
    base_path = os.path.join(config['output']['output_path'], "Target-binding-metrics")
    # Check if the directory exists
    if not os.path.exists(base_path):
        return None

    def calculate_mean_from_file(file_path):
        """
        Extract the values from the second column of the file and calculate the average, ignoring the header.
        :param file_path: path
        :return: mean value or 'NA'
        """
        if not os.path.exists(file_path):
            return 'NA'

        try:
            data = pd.read_csv(file_path, sep='\t', header=0, usecols=[1])
            if data.empty:
                return 'NA'

            mean_value = data.iloc[:, 0].mean()
            return f"{mean_value:.9f}"
        except Exception as e:
            #print(f"Error processing file {file_path}: {e}")
            return 'NA'

    def determine_file_prefix():
        """
        Determine whether files in the directory start with 'glide' or 'vina'
        :return: file name prefix ('glide' or 'vina')
        """
        files = os.listdir(base_path)
        #print(f"Files in the directory: {files}")  # Debugging information
        glide_files = [
            'glide-docking-score-top-10.txt',
            'glide-fraction-of-PM.txt',
            'glide-in-situ-docking-score-PM.txt',
            'glide-in-situ-docking-score-top-10.txt',
            'glide-openeye-rmsd-score.txt',
            'glide-rdkit-rmsd-score.txt',
            'glide-openeye-rmsd-score-top10.txt',
            'glide-rdkit-rmsd-score-top10.txt'
        ]
        vina_files = [
            'vina-docking-score-top-10.txt',
            'vina-fraction-of-PM.txt',
            'vina-in-situ-docking-score-PM.txt',
            'vina-in-situ-docking-score-top-10.txt',
            'vina-openeye-rmsd-score.txt',
            'vina-rdkit-rmsd-score.txt',
            'vina-openeye-rmsd-score-top10.txt',
            'vina-rdkit-rmsd-score-top10.txt'
        ]

        for file in files:
            if file in glide_files:
                #print("Detected prefix: glide")  # Debugging information
                return 'glide'
            elif file in vina_files:
                #print("Detected prefix: vina")  # Debugging information
                return 'vina'
        return None

    prefix = determine_file_prefix()
    if prefix is None:
        #print("No glide or vina files found in the directory.")
        return

    file_names = [
        f'{prefix}-docking-score-top-10.txt',
        f'{prefix}-fraction-of-PM.txt',
        f'{prefix}-in-situ-docking-score-PM.txt',
        f'{prefix}-in-situ-docking-score-top-10.txt',
        f'{prefix}-openeye-rmsd-score.txt',
        f'{prefix}-rdkit-rmsd-score.txt',
        f'{prefix}-openeye-rmsd-score-top10.txt',
        f'{prefix}-rdkit-rmsd-score-top10.txt'
    ]

    results = {}
    for file_name in file_names:
        file_path = os.path.join(base_path, file_name)
        if os.path.exists(file_path):
            mean_value = calculate_mean_from_file(file_path)
            results[file_name.replace('.txt', '')] = mean_value
        else:
            results[file_name.replace('.txt', '')] = 'NA'

    output_file_path = os.path.join(base_path, 'Target-binding-metrics.txt')
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    with open(output_file_path, 'w') as output_file:
        output_file.write("Target-binding-metrics\tMean Values\n")
        for key, value in results.items():
            if value != 'NA':  # Only write detected files
                output_file.write(f"{key}\t{value}\n")
