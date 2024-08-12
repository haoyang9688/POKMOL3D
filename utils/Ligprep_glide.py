import os
import subprocess
from multiprocessing import Pool

def process_sdf_file(input_sdf_file, output_sdf_file):
    ligprep_command = f"ligprep -isd {input_sdf_file} -epik -ph 7.0 -pht 0.2 -s 1 -NJOBS 12 -omae {output_sdf_file}.maegz"
    result = subprocess.run(ligprep_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        print(f" mae: {output_sdf_file}.maegz")
    else:
        print(f"mae error: {input_sdf_file}")
        #print("error log：")
        #print(result.stderr)

def process_folder(folder_path, output_base_path):
    sdf_files = [f for f in os.listdir(folder_path) if f.endswith(".sdf") and os.path.isfile(os.path.join(folder_path, f))]
    if sdf_files:
        for sdf_file in sdf_files:
            input_sdf_file = os.path.join(folder_path, sdf_file)
            folder_index = int(os.path.basename(folder_path)) % 32
            output_folder = os.path.join(output_base_path, str(folder_index), 'Ligprep-output-maegz')
            if not os.path.exists(output_folder):
                os.makedirs(output_folder, exist_ok=True)  
            output_sdf_file = os.path.join(output_folder, os.path.splitext(sdf_file)[0])
            process_sdf_file(input_sdf_file, output_sdf_file)

def glide_ligprep_process(config):
    config_output_path = config['output']['output_path']
    main_directory = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
    ligprep_output_path = config['settings']['Redocking_settings']['prepared_ligands_path']
    os.makedirs(ligprep_output_path, exist_ok=True)

    folders = [os.path.join(main_directory, str(i)) for i in range(32) if os.path.isdir(os.path.join(main_directory, str(i)))]

    with Pool(processes=100) as pool:
        pool.starmap(process_folder, [(folder, ligprep_output_path) for folder in folders])

