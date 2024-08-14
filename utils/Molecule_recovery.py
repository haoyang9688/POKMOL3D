import os
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def calculate_similarity(std_smiles, gen_smiles):
    try:
        std_mol = Chem.MolFromSmiles(std_smiles)
        gen_mol = Chem.MolFromSmiles(gen_smiles)
        std_fp = AllChem.GetMorganFingerprintAsBitVect(std_mol, 2, nBits=2048)
        gen_fp = AllChem.GetMorganFingerprintAsBitVect(gen_mol, 2, nBits=2048)
        score = DataStructs.TanimotoSimilarity(std_fp, gen_fp)
        return score
    except Exception as e:
        print(f"Error calculating similarity: {e}")
        return 0.0

def calculate_molecule_recovery(config):
    active_data_path = config['active']['active_data_path']
    base_path = config['output']['output_path']
    thresholds = config['settings']['Active_recovery_settings']['similarity_thresholds']

    output_base_folder = os.path.join(base_path, 'Recovery-active-metrics', 'Molecule')
    os.makedirs(output_base_folder, exist_ok=True)

    for folder_number in range(32):
        gen_folder_path = os.path.join(base_path, 'General-molecular-quality-metrics', 'Uniqueness', str(folder_number))
        std_folder_path = os.path.join(active_data_path, str(folder_number))
        
        if not os.path.exists(gen_folder_path):
            #print(f": {gen_folder_path}")
            continue
        if not os.path.exists(std_folder_path):
            #print(f"The standard file path does not exist: {std_folder_path}")
            continue
        
        gen_files = [f for f in os.listdir(gen_folder_path) if f.endswith('.txt')]
        std_files = [f for f in os.listdir(std_folder_path) if f.endswith('.txt')]

        # Create subfolders 0-31 in the output folder
        output_folder_path = os.path.join(output_base_folder, str(folder_number))
        os.makedirs(output_folder_path, exist_ok=True)
        
        all_similarity_scores_files = {threshold: open(os.path.join(output_folder_path, f'Molecule-recovery-{threshold}.txt'), 'w') for threshold in thresholds}
        
        for gen_file_name in gen_files:
            gen_file_path = os.path.join(gen_folder_path, gen_file_name)
            with open(gen_file_path, 'r') as gen_file:
                gen_smiles_lines = gen_file.readlines()
                
                for std_file_name in std_files:
                    std_txt_file_path = os.path.join(std_folder_path, std_file_name)
                    with open(std_txt_file_path, 'r') as std_file:
                        for std_smiles_line in std_file:
                            std_id, std_smiles = std_smiles_line.strip().split('\t', 1)
                            for gen_smiles_line in gen_smiles_lines:
                                gen_id, gen_smiles = gen_smiles_line.strip().split('\t', 1)
                                similarity_score = calculate_similarity(std_smiles, gen_smiles)
                                for threshold in thresholds:
                                    if similarity_score >= threshold:
                                        all_similarity_scores_files[threshold].write(f"{std_id}\t{gen_id}\t{std_smiles}\t{gen_smiles}\t{similarity_score}\n")
                                        break

        for file in all_similarity_scores_files.values():
            file.close()

        #print(f"molecule-recovery {folder_number}")

def load_config(config_path='config.yml'):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

# 
if __name__ == "__main__":
    config = load_config('config.yml')
    if config['evaluate']['active_recovery']:
        calculate_molecule_recovery(config)
