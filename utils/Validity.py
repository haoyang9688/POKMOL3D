import os
import re
from rdkit import Chem

def calculate_validity(config):
    
    data_path = config['data']['input_path']
    output_path = config['output']['output_path']

    
    for folder_number in range(32):
        folder_path = os.path.join(data_path, str(folder_number))

       
        if not os.path.isdir(folder_path):
            continue 

        folder_number = os.path.basename(folder_path)

        general_molecular_quality_folder = os.path.join(output_path, 'General-molecular-quality-metrics')
        os.makedirs(general_molecular_quality_folder, exist_ok=True)

        validity_folder_path = os.path.join(general_molecular_quality_folder, 'Validity')
        os.makedirs(validity_folder_path, exist_ok=True)

        subfolder_path = os.path.join(validity_folder_path, folder_number)
        os.makedirs(subfolder_path, exist_ok=True)

        sdf_files = [f for f in os.listdir(folder_path) if f.endswith('.sdf')]
        smiles_data = []

        for sdf_file in sdf_files:
            sdf_file_path = os.path.join(folder_path, sdf_file)

            
            match = re.match(r'([0-9]+)\.sdf', sdf_file)
            if match:
                file_id = int(match.group(1))  
            else:
                continue

           
            mol_supplier = Chem.SDMolSupplier(sdf_file_path)
            if mol_supplier is not None:
    
                for idx, mol in enumerate(mol_supplier):
                    if mol is not None:
                        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                        smiles_data.append((file_id, canonical_smiles))
                    else:
                        print(f"{sdf_file} {idx + 1} error")
            else:
                print(f"{sdf_file} error")

        smiles_data.sort()
        output_file_path = os.path.join(subfolder_path, 'validity.txt')
        with open(output_file_path, 'w') as output_file:
            for file_id, canonical_smiles in smiles_data:
                output_file.write(f"{file_id}\t{canonical_smiles}\n")

