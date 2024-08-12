import os
import subprocess
from multiprocessing import Pool

def process_sdf_file(input_sdf_file, output_sdf_file):
    ligprep_command = f"ligprep -isd {input_sdf_file} -bff 16 -g -i 0 -nt -s 1 -NJOBS 12 -osd {output_sdf_file}"
    subprocess.run(ligprep_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

def process_folder(folder_path, output_sdf_path):
    sdf_files = [f for f in os.listdir(folder_path) if f.endswith(".sdf") and os.path.isfile(os.path.join(folder_path, f))]
    for sdf_file in sdf_files:
        input_sdf_file = os.path.join(folder_path, sdf_file)
        output_folder = os.path.join(output_sdf_path, os.path.basename(folder_path))
        os.makedirs(output_folder, exist_ok=True)
        output_sdf_file = os.path.join(output_folder, sdf_file)
        process_sdf_file(input_sdf_file, output_sdf_file)

def calculate_OPLS3_preparation(config):
    config_output_path = config['output']['output_path']
    config_input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
    output_sdf_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Force-field-molecules')
    os.makedirs(output_sdf_path, exist_ok=True)
    folders = [os.path.join(config_input_path, str(i)) for i in range(32) if os.path.exists(os.path.join(config_input_path, str(i)))]
    
    with Pool(processes=120) as pool:
        pool.starmap(process_folder, [(folder, output_sdf_path) for folder in folders])
    