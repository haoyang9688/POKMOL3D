import os
import subprocess
from multiprocessing import Pool

#Define a function to check if there is a score in the result file
def check_for_scores(output_file):
    with open(output_file, 'r') as result_file:
        for line in result_file:
            if line.startswith("   1 "):
                return True
    return False

#Define a function to extract the highest score from the result file
def extract_best_score(output_file):
    best_score = float('-inf')
    with open(output_file, 'r') as result_file:
        lines = result_file.readlines()
        for line in lines:
            if line.startswith("   1 "):
                score = float(line.split()[1])
                best_score = max(best_score, score)
    return best_score

def process_task(task):
    numeric_folder, receptor_file, ligand_file, docking_recepter, docking_ligand, output_base_folder, qvina_executable = task

    numeric_folder_path = os.path.join(docking_recepter, numeric_folder)
    receptor_pdbqt_file = os.path.join(numeric_folder_path, receptor_file)

    output_folder = os.path.join(output_base_folder, numeric_folder)
    output_pdbqt_folder = os.path.join(output_folder)
    qvina_text_folder = os.path.join(output_folder, "qvina-text")

    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(output_pdbqt_folder, exist_ok=True)
    os.makedirs(qvina_text_folder, exist_ok=True)

    ligand_folder_path = os.path.join(docking_ligand, numeric_folder, "Ligprep-output-pdbqt")
    ligand_pdbqt_file = os.path.join(ligand_folder_path, ligand_file)

    output_pdbqt = os.path.join(output_pdbqt_folder, f"output-{ligand_file}")
    qvina_output_file = os.path.join(qvina_text_folder, f"qvina-out-{ligand_file}.txt")
    config_file = os.path.join(numeric_folder_path, "docking-config.txt")

    try:
        qvina_command = [qvina_executable, '--config', config_file, '--score_only', '--receptor', receptor_pdbqt_file, '--ligand', ligand_pdbqt_file, '--out', output_pdbqt , '--cpu', '1']
        with open(qvina_output_file, 'w') as output_file:
            subprocess.run(qvina_command, stdout=output_file, stderr=subprocess.STDOUT)
        
        if not check_for_scores(qvina_output_file):
            #print(f"No score, give up docking: {ligand_file}")
            return
    
    except Exception as e:
        #print(f"error: {e}")
        pass
def vina_in_situ_docking_main(config):
    docking_recepter = os.path.join(config['POKMOL3D_path'], 'Source/vina-prepared-receptors')
    config_output_path = config['output']['output_path']
    docking_ligand = os.path.join(config_output_path, "Target-binding-metrics", 'Prepared-add-h-ligands')
    output_base_folder = os.path.join(config_output_path, 'Target-binding-metrics', 'Vina-In-situ-Docking-results')
    os.makedirs(output_base_folder, exist_ok=True)
    qvina_executable_path = config['docking_env']['vina_path']
    qvina_executable = os.path.join(qvina_executable_path, "qvina2")

    tasks = []
    for numeric_folder in os.listdir(docking_recepter):
        numeric_folder_path = os.path.join(docking_recepter, numeric_folder)
        if os.path.isdir(numeric_folder_path):
            receptor_files = [f for f in os.listdir(numeric_folder_path) if f.endswith(".pdbqt")]
            ligand_folder_path = os.path.join(docking_ligand, numeric_folder, "Ligprep-output-pdbqt")
            if os.path.exists(ligand_folder_path):
                ligand_files = [f for f in os.listdir(ligand_folder_path) if f.endswith(".pdbqt")]
                for receptor_file in receptor_files:
                    for ligand_file in ligand_files:
                        tasks.append((numeric_folder, receptor_file, ligand_file, docking_recepter, docking_ligand, output_base_folder, qvina_executable))

    with Pool(processes=100) as pool:
        pool.map(process_task, tasks)
