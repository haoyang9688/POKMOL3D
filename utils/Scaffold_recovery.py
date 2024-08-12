import os
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, DataStructs
import yaml

def load_config(config_file='config.yml'):
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def process_recovery_of_active_molecules(std_smiles, gen_smiles):
    try:
        std_mol = Chem.MolFromSmiles(std_smiles)
        gen_mol = Chem.MolFromSmiles(gen_smiles)
        
        std_scaffold = MurckoScaffold.GetScaffoldForMol(std_mol)
        gen_scaffold = MurckoScaffold.GetScaffoldForMol(gen_mol)
        
        std_scaffold_mol = Chem.MolFromSmiles(Chem.MolToSmiles(std_scaffold))
        gen_scaffold_mol = Chem.MolFromSmiles(Chem.MolToSmiles(gen_scaffold))
        
        std_fp = AllChem.GetMorganFingerprintAsBitVect(std_scaffold_mol, radius=2, nBits=2048)
        gen_fp = AllChem.GetMorganFingerprintAsBitVect(gen_scaffold_mol, radius=2, nBits=2048)
        
        score = DataStructs.TanimotoSimilarity(std_fp, gen_fp)

        return score
    except Exception as e:
        print(f"Error processing molecules: {e}")
        return 0.0

def calculate_scafold_recovery(config):
    active_data_path = config['active']['active_data_path']
    base_path = config['output']['output_path']
    thresholds = config['settings']['Active_recovery_settings']['similarity_thresholds']
    
    #Get output path
    scaffold_recovery_output_base = os.path.join(base_path, 'Recovery-active-metrics', 'Scaffold')
    os.makedirs(scaffold_recovery_output_base, exist_ok=True)
    
    for folder_number in range(32):
        gen_folder_path = os.path.join(base_path, 'General-molecular-quality-metrics', 'Uniqueness', str(folder_number))
        std_folder_path = os.path.join(active_data_path, str(folder_number))
        
        if not os.path.exists(gen_folder_path):
            #print(f"The generated file path does not exist: {gen_folder_path}")
            continue
        if not os.path.exists(std_folder_path):
            #print(f"The generated file path does not exist: {std_folder_path}")
            continue
        
        gen_files = [f for f in os.listdir(gen_folder_path) if f.endswith('.txt')]
        std_files = [f for f in os.listdir(std_folder_path) if f.endswith('.txt')]

        #Create subfolders 0-31 in the Recovery_scaffolds metrics folder
        output_folder_path = os.path.join(scaffold_recovery_output_base, str(folder_number))
        os.makedirs(output_folder_path, exist_ok=True)
        
        for threshold in thresholds:
            all_similarity_scores_file_path = os.path.join(output_folder_path, f'Scaffold-recovery-{threshold}.txt')
            with open(all_similarity_scores_file_path, 'w') as f_all_similarity_scores:
                for std_file_name in std_files:
                    std_txt_file_path = os.path.join(std_folder_path, std_file_name)
                    with open(std_txt_file_path, 'r') as std_file:
                        std_data = std_file.read().splitlines()
                        for std_line in std_data:
                            std_id, std_smiles = std_line.split('\t')
                            for gen_file_name in gen_files:
                                gen_txt_file_path = os.path.join(gen_folder_path, gen_file_name)
                                with open(gen_txt_file_path, 'r') as gen_file:
                                    gen_smiles_data = gen_file.read().splitlines()
                                    for gen_smiles_line in gen_smiles_data:
                                        gen_mol_id, gen_smiles = gen_smiles_line.split('\t')
                                        murcko_similarity_score = process_recovery_of_active_molecules(std_smiles, gen_smiles)
                                        if murcko_similarity_score >= threshold:
                                            f_all_similarity_scores.write(f"{std_id}\t{gen_mol_id}\t{std_smiles}\t{gen_smiles}\t{murcko_similarity_score}\n")
                                            break  #End the current loop and continue with the next standard molecule
            

# 主程序入口
if __name__ == "__main__":
    config = load_config('config.yml')
    calculate_scafold_recovery(config)
