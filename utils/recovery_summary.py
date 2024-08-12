import os

def count_first_column_lines(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        return len([line for line in lines if line.strip().split('\t')[0]])

def process_folders(base_dir_1, base_dir_2, output_file):
    results = []
    ratios = []

    for i in range(32):  
        folder_1 = os.path.join(base_dir_1, str(i))
        folder_2 = os.path.join(base_dir_2, str(i))

        # Check if the folder exists, skip if it does not exist
        if not os.path.exists(folder_1) or not os.path.exists(folder_2):
            continue

        total_lines_1 = 0
        total_lines_2 = 0

        for file_name in os.listdir(folder_1):
            if file_name.endswith('.txt'):
                file_path_1 = os.path.join(folder_1, file_name)
                lines_1 = count_first_column_lines(file_path_1)
                total_lines_1 += lines_1

        for file_name in os.listdir(folder_2):
            if file_name.endswith('.txt'):
                file_path_2 = os.path.join(folder_2, file_name)
                lines_2 = count_first_column_lines(file_path_2)
                total_lines_2 += lines_2

        if total_lines_2 != 0:
            ratio = total_lines_1 / total_lines_2
        else:
            ratio = float('inf') 

        results.append(f"Target\t{i}\trecovery:\t{ratio:.3f}\n")
        ratios.append(ratio)

    mean_ratio = sum(ratios) / len(ratios) if ratios else 0

    results.append(f"mean:\t{mean_ratio:.3f}\n")

    with open(output_file, 'w') as out_file:
        out_file.writelines(results)

def extract_mean_recoveries(molecule_file, scaffold_file, output_file):
    with open(molecule_file, 'r') as mol_file:
        last_line_mol = mol_file.readlines()[-1].strip()
        mol_mean_recovery = last_line_mol.split(':')[-1].strip()

    with open(scaffold_file, 'r') as scaf_file:
        last_line_scaf = scaf_file.readlines()[-1].strip()
        scaf_mean_recovery = last_line_scaf.split(':')[-1].strip()

    with open(output_file, 'w') as out_file:
        out_file.write(f"Molecule mean recovery:\t{mol_mean_recovery}\n")
        out_file.write(f"Scaffold mean recovery:\t{scaf_mean_recovery}\n")

def calculate_recovery(config):
    output_base_path = config['output']['output_path']
    base_dir_2 = config['active']['active_data_path']

    base_dir_1_molecule = os.path.join(output_base_path, "Recovery-active-metrics", "Molecule")
    output_file_molecule = os.path.join(base_dir_1_molecule, "Molecule-recovery.txt")

    base_dir_1_scaffold = os.path.join(output_base_path, "Recovery-active-metrics", "Scaffold")
    output_file_scaffold = os.path.join(base_dir_1_scaffold, "Scaffold-recovery.txt")

    process_folders(base_dir_1_molecule, base_dir_2, output_file_molecule)
    process_folders(base_dir_1_scaffold, base_dir_2, output_file_scaffold)

    # Extract and write mean recoveries to a new file
    active_recovery_metrics_file = os.path.join(output_base_path, "Recovery-active-metrics", "Active-recovery-metrics.txt")
    extract_mean_recoveries(output_file_molecule, output_file_scaffold, active_recovery_metrics_file)
