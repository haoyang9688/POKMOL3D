import os
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from itertools import combinations
from rdkit.Chem import AllChem
from multiprocessing import Pool

def calculate_scaffold_similarity(gen_smiles):
    """Extract the basic skeleton of molecules and calculate fingerprints"""
    gen_mol = Chem.MolFromSmiles(gen_smiles)
    if gen_mol is None or gen_mol.GetRingInfo().NumRings() == 0 or gen_mol.GetNumAtoms() < 3:
        return None

    gen_scaffold = MurckoScaffold.GetScaffoldForMol(gen_mol)
    if gen_scaffold is None:
        return None

    gen_scaffold_mol = Chem.MolFromSmiles(Chem.MolToSmiles(gen_scaffold))
    gen_fp = AllChem.GetMorganFingerprintAsBitVect(gen_scaffold_mol, radius=2, nBits=2048)

    return gen_fp

def process_single_folder(config, numeric_subfolder):
    config_output_path = config['output']['output_path']
    input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics')
    uniqueness_base_folder = os.path.join(input_path, 'Uniqueness')
    diversity_base_folder = os.path.join(input_path, 'Scaffold-Interal-diversity')
    os.makedirs(diversity_base_folder, exist_ok=True)

    numeric_subfolder_path = os.path.join(uniqueness_base_folder, str(numeric_subfolder))
    if not os.path.isdir(numeric_subfolder_path):
        return

    smiles_list = []
    file_paths = [os.path.join(numeric_subfolder_path, file_name) for file_name in os.listdir(numeric_subfolder_path) if file_name.endswith('.txt')]

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                smiles = parts[1].strip()
                smiles_list.append(smiles)

    output_folder = os.path.join(diversity_base_folder, str(numeric_subfolder))
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, 'scaffold-interal-diversity.txt')

    try:
        with open(output_path, 'w') as output_file:
            tanimoto_sum = 0.0
            non_empty_molecules = 0

            for i, j in combinations(range(len(smiles_list)), 2):
                gen_smiles_i, gen_smiles_j = smiles_list[i], smiles_list[j]
                gen_fp_i = calculate_scaffold_similarity(gen_smiles_i)
                gen_fp_j = calculate_scaffold_similarity(gen_smiles_j)
                if gen_fp_i is None or gen_fp_j is None:
                    continue
                similarity_score = DataStructs.TanimotoSimilarity(gen_fp_i, gen_fp_j)
                tanimoto_sum += similarity_score
                non_empty_molecules += 1

        total_molecules = len(smiles_list)
        overall_diversity = 1 - (tanimoto_sum / total_molecules**2) if total_molecules > 0 else 0.0

        with open(output_path, 'a') as overall_output_file:
            overall_output_file.write(f"{os.path.basename(numeric_subfolder_path)}\tvalue:\t{overall_diversity}\n")

    except IOError as e:
        print(f"Error: {e}")

def calculate_Scaffold_Tdiv(config):
    subfolder_range = range(0, 32)
    with Pool() as pool:
        pool.starmap(process_single_folder, [(config, i) for i in subfolder_range])
