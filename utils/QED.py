import os
import glob
from rdkit import Chem
from rdkit.Chem import QED
from multiprocessing import Pool


def process_txt_file(txt_file):
    with open(txt_file, 'r') as f:
        lines = f.readlines()

    results = []

    for line in lines:
        data = line.strip().split('\t')
        if len(data) < 2:
            #print(f"Warning: Skip lines and files with insufficient data:{txt_file}")
            continue

        molecule_id, smiles = data[0], data[1]
        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            qed_value = QED.qed(mol)
            results.append((molecule_id, smiles, qed_value))
        else:
            print(f"Warning: Skip lines and files with insufficient data{txt_file}")

    return results

def calculate_QED(config):
    config_output_path = config['output']['output_path']
    input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics')
    uniqueness_base_folder = os.path.join(input_path, 'Uniqueness')
    output_base_folder = os.path.join(input_path, 'QED')
    os.makedirs(output_base_folder, exist_ok=True)
    
    #print(f"Output base folder: {output_base_folder}")  

    for i in range(32):
        numeric_subfolder = str(i)
        uniqueness_folder = os.path.join(uniqueness_base_folder, numeric_subfolder)
        output_folder = os.path.join(output_base_folder, numeric_subfolder)
        if not os.path.exists(uniqueness_folder):
            #print(f" {uniqueness_folder} non-existent. Skip.")
            continue
        os.makedirs(output_folder, exist_ok=True)

        txt_files = glob.glob(os.path.join(uniqueness_folder, '*.txt'))

        if not txt_files:
            #print(f" {uniqueness_folder} not contain any. txt files. skip。")
            continue

        #print(f"handel {numeric_subfolder} 'Uniqueness' ...")

        summary_file = os.path.join(output_folder, 'QED.txt')
        with open(summary_file, 'w') as summary_f:
            summary_f.write("ID\tSMILES\tQED\n")

        with Pool() as pool:
            results = pool.map(process_txt_file, txt_files)
            flattened_results = [result for sublist in results for result in sublist]  
            with open(summary_file, 'a') as summary_f:
                for result in flattened_results:
                    summary_f.write('\t'.join(map(str, result)) + '\n')



