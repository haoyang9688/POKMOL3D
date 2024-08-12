import os
from rdkit import Chem

def calculate_usability(config):
    try:
        #Effective set of elements
        valid_element_set = {'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'H'}
        config_output_path = config['output']['output_path']
        #Root directory path
        root_folder = os.path.join(config_output_path, 'General-molecular-quality-metrics')
        
        #Unique folder path
        uniqueness_folder = os.path.join(root_folder, 'Uniqueness')
    
        #Create Usability folder
        usability_folder = os.path.join(root_folder, 'Usability')
        os.makedirs(usability_folder, exist_ok=True)
        
        #Traverse folders 0-31
        for i in range(32):
            subfolder_path = os.path.join(uniqueness_folder, str(i))
            unique_file = os.path.join(subfolder_path, 'unique.txt')
            
            if not os.path.isfile(unique_file):
                continue
            
            #Create an empty list to store all standardized SMILES strings and their corresponding ID numbers in the current subfolders
            all_valid_smiles = []
            
            #Read the contents of the unique.txt file
            with open(unique_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        id_number, smiles_string = parts
                        #Verify the validity of SMILES strings using RDKit
                        mol = Chem.MolFromSmiles(smiles_string, sanitize=False)
                        if mol is not None:
                            #Extract element symbols contained in molecules
                            elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
                            unwanted_elements = elements - valid_element_set
                            if not unwanted_elements:
                                all_valid_smiles.append((id_number, smiles_string))
            
            sub_usability_folder = os.path.join(usability_folder, str(i))
            os.makedirs(sub_usability_folder, exist_ok=True)
            
            output_file = os.path.join(sub_usability_folder, 'usability.txt')
            with open(output_file, 'w') as f:
                for id_number, valid_smi in all_valid_smiles:
                    f.write(f"{id_number}\t{valid_smi}\n")
            
            #print(f"Valid SMILES have been saved to {output_file}")
    
    except Exception as e:
        print(f"Usability Calculation failed: {e}")


