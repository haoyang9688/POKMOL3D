import os
import subprocess
import glob
from multiprocessing import Pool

def process_folder(args):
    return calculate_redocking(*args)

def calculate_redocking(number_folder, docking_recepter_cwd, docking_ligand_cwd, output_base_path, schrodinger_path):
    #Build the path of the current digital folder
    number_folder_path = os.path.join(docking_recepter_cwd, number_folder)
    #Add or update receptor pathways
    receptor_folder = os.path.join(number_folder_path, "glide-grid.zip")

    #Create or retrieve output directory
    target_binding_metrics_folder = os.path.join(output_base_path, "Target-binding-metrics")
    docking_results_folder = os.path.join(target_binding_metrics_folder, "Glide-Docking-results", number_folder)
    if not os.path.exists(docking_results_folder):
        os.makedirs(target_binding_metrics_folder, exist_ok=True)

    #Retrieve the paths of all. in files in the current digital folder
    docking_config_paths = glob.glob(os.path.join(number_folder_path, "*.in"))
    if not docking_config_paths:
        print(f"No .in file found in {number_folder_path}")
        return
    docking_config_path = docking_config_paths[0]

    #Retrieve the paths of all ligand files in the current digital folder
    ligand_folder = os.path.join(docking_ligand_cwd, number_folder, "Ligprep-output-maegz")
    if not os.path.exists(ligand_folder):
        print(f"Ligand folder does not exist: {ligand_folder}")
        return
    ligand_files = [f for f in os.listdir(ligand_folder) if f.endswith(".maegz")]
    for ligand_file in ligand_files:
        ligand_path = os.path.join(ligand_folder, ligand_file)
        prefix = ligand_file.split(".maegz")[0]
        new_docking_config_name = f"{ligand_file.split('.')[0]}.in"
        new_docking_config_path = os.path.join(number_folder_path, new_docking_config_name)

        #Rename the original. in file
        os.rename(docking_config_path, new_docking_config_path)
        docking_config_path = new_docking_config_path

        #Modify parameters in the. in file
        with open(new_docking_config_path, 'r') as f:
            lines = f.readlines()

        modified_lines = []
        for line in lines:
            if line.startswith("GRIDFILE"):
                line = f"GRIDFILE     {receptor_folder}\n"
            elif line.startswith("LIGANDFILE"):
                line = f"LIGANDFILE   {ligand_path}\n"
            elif line.startswith("OUTPUTDIR"):
                line = f"OUTPUTDIR    {docking_results_folder}\n"
            modified_lines.append(line)

        #Write the modified lines back to the. in file
        with open(new_docking_config_path, 'w') as f:
            f.writelines(modified_lines)

       #Run Glide command
        glide_executable = os.path.join(schrodinger_path, "glide")
        docking_command = [glide_executable, new_docking_config_name]
        result = subprocess.run(docking_command, cwd=number_folder_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

def glide_redocking_main(config):
    docking_recepter_cwd = os.path.join(config['POKMOL3D_path'], 'Source/glide-prepared-receptors')
    docking_ligand_cwd = config['settings']['Redocking_settings']['prepared_ligands_path']
    config_output_path = config['output']['output_path']
    schrodinger_path = config['docking_env']['Schrodinger_path']
    output_base_path = os.path.join(config_output_path)
    number_folders = [f for f in os.listdir(docking_ligand_cwd) if os.path.isdir(os.path.join(docking_ligand_cwd, f)) and f.isdigit()]
    args = [(folder, docking_recepter_cwd, docking_ligand_cwd, output_base_path, schrodinger_path) for folder in number_folders]

    with Pool(processes=100) as pool:
        pool.map(process_folder, args)


