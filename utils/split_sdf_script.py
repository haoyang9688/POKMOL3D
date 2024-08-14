import os
from rdkit import Chem

def split_and_save_sdf(input_sdf_path, output_folder):
    suppl = Chem.SDMolSupplier(input_sdf_path)
    os.makedirs(output_folder, exist_ok=True)
    
    for idx, mol in enumerate(suppl, start=1):
        if mol is not None:
            output_sdf_path = os.path.join(output_folder, f'{idx}.sdf')
            with open(output_sdf_path, 'w') as sdf_file:
                sdf_file.write(Chem.MolToMolBlock(mol))

if __name__ == "__main__":
    root_directory = '/your/Model/path'
    all_subfolders = [f for f in os.listdir(root_directory) if os.path.isdir(os.path.join(root_directory, f))]

    for subfolder in all_subfolders:
        subfolder_path = os.path.join(root_directory, subfolder)
        numeric_subfolders = [f for f in os.listdir(subfolder_path) if f.isdigit()]
        for numeric_subfolder in numeric_subfolders:
            numeric_subfolder_path = os.path.join(subfolder_path, numeric_subfolder)
            sdf_files = [f for f in os.listdir(numeric_subfolder_path) if f.endswith('.sdf')]
            sdf_folder_path = os.path.join(numeric_subfolder_path, 'Split-SDF')
            os.makedirs(sdf_folder_path, exist_ok=True)
            for sdf_file in sdf_files:
                input_sdf_path = os.path.join(numeric_subfolder_path, sdf_file)
                output_folder = sdf_folder_path
                split_and_save_sdf(input_sdf_path, output_folder)
