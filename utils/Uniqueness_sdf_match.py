import os
import shutil
from rdkit import Chem

def match_and_copy_sdf_files(config):
    config_output_path = config['output']['output_path']
    base_txt_path = os.path.join(config_output_path, "General-molecular-quality-metrics", "Uniqueness")
    base_sdf_path = config['data']['input_path']
    output_base_path = os.path.join(config_output_path, "General-molecular-quality-metrics", "Uniqueness-sdf")

    for folder_number in range(32):
        txt_folder_path = os.path.join(base_txt_path, str(folder_number))
        sdf_folder_path = os.path.join(base_sdf_path, str(folder_number))

        # Only proceed if both the txt_folder_path and sdf_folder_path exist
        if os.path.isdir(txt_folder_path) and os.path.isdir(sdf_folder_path):
            output_folder_path = os.path.join(output_base_path, str(folder_number))
            os.makedirs(output_folder_path, exist_ok=True)

            # Process each .txt file in the txt_folder_path
            for txt_file_name in os.listdir(txt_folder_path):
                if txt_file_name.endswith('.txt'):
                    txt_file_path = os.path.join(txt_folder_path, txt_file_name)
                    with open(txt_file_path, 'r') as txt_file:
                        for line in txt_file:
                            parts = line.strip().split()
                            if len(parts) >= 1:
                                sdf_index = parts[0]
                                sdf_file_name = f"{sdf_index}.sdf"
                                sdf_file_path = os.path.join(sdf_folder_path, sdf_file_name)
                                
                                if os.path.exists(sdf_file_path):
                                    shutil.copy(sdf_file_path, output_folder_path)
                                else:
                                    print(f"SDF file {sdf_file_path} does not exist.")
        else:
            print(f"Folder {txt_folder_path} or {sdf_folder_path} does not exist.")

