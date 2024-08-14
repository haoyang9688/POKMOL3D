import os
import gzip
import shutil
from multiprocessing import Pool
import yaml

def process_folder(number_folder, res_omaegz_root):
    number_folder_path = os.path.join(res_omaegz_root, number_folder)
    docking_txt_folder = os.path.join(res_omaegz_root, number_folder)

    if os.path.exists(number_folder_path) and os.path.isdir(number_folder_path):
    
        docking_txt_path = os.path.join(docking_txt_folder, "in-situ-docking-scores.txt")
        if not os.path.exists(docking_txt_folder):
            os.makedirs(docking_txt_folder)
        
        for sdfgz_file_name in os.listdir(number_folder_path):
            if sdfgz_file_name.endswith(".sdfgz"):
                sdfgz_file_path = os.path.join(number_folder_path, sdfgz_file_name)

                with gzip.open(sdfgz_file_path, 'rb') as f_in:
                    sdf_file_name = os.path.splitext(sdfgz_file_name)[0] + ".sdf"
                    sdf_file_path = os.path.join(number_folder_path, sdf_file_name)
                    with open(sdf_file_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                
                with open(sdf_file_path, 'r') as sdf_file:
                    lines = sdf_file.readlines()
                    docking_score = None
                    for i in range(len(lines)):
                        if lines[i].strip() == "> <r_i_docking_score>":
                            docking_score = lines[i+1].strip()
                            break
                
                if docking_score:
                    sdf_prefix = os.path.splitext(sdfgz_file_name)[0].split('_')[0]  
                    with open(docking_txt_path, 'a') as docking_txt_file:
                        docking_txt_file.write(f"{sdf_prefix}: {docking_score}\n")

def uzip_sdfgz_main(config):
    config_output_path = config['output']['output_path']
    res_omaegz_root = os.path.join(config_output_path, "Target-binding-metrics", "In-situ-Docking-results")
    
    number_folders = [folder for folder in os.listdir(res_omaegz_root) if folder.isdigit() and 0 <= int(folder) <= 31]

    with Pool() as pool:
        pool.starmap(process_folder, [(folder, res_omaegz_root) for folder in number_folders])

if __name__ == "__main__":
    try:
    
        with open("config.yaml", "r") as file:
            config = yaml.safe_load(file)
        
        uzip_sdfgz_main(config)
    except Exception as e:
       pass
