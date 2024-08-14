import os
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

#Predefine all possible key types
all_bond_keys = ['C-C', 'C=C', 'C≡C', 'C-N', 'C=N', 'C≡N', 'C-O', 'C=O', 'C:O']

#Calculate bond length information in molecules
def calculate_bond_lengths(mol):
    conf = mol.GetConformer()
    bond_lengths = {key: [] for key in all_bond_keys}
    
    for bond in mol.GetBonds():
        atom_idx_1 = bond.GetBeginAtomIdx()
        atom_idx_2 = bond.GetEndAtomIdx()
        bond_length = rdMolTransforms.GetBondLength(conf, atom_idx_1, atom_idx_2)

        atom1_type = mol.GetAtomWithIdx(atom_idx_1).GetSymbol()
        atom2_type = mol.GetAtomWithIdx(atom_idx_2).GetSymbol()

        if atom1_type > atom2_type:
            atom1_type, atom2_type = atom2_type, atom1_type
        
        bond_key = f'{atom1_type}-{atom2_type}'
        
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            pass
        elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            bond_key = f'{atom1_type}={atom2_type}'
        elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            bond_key = f'{atom1_type}≡{atom2_type}'
        elif bond.GetBeginAtom().GetIsAromatic() and bond.GetEndAtom().GetIsAromatic():
            bond_key = f'{atom1_type}:{atom2_type}'
        
        if bond_key in bond_lengths:
            bond_lengths[bond_key].append(bond_length)
    
    return bond_lengths

def process(data_path_1, data_path_2, folder_name, output_path_1, output_path_2):
    #No optimized molecular data path
    folder_path_1 = os.path.join(data_path_1, folder_name)
    bond_length_folder_1 = os.path.join(output_path_1, folder_name)
    os.makedirs(bond_length_folder_1, exist_ok=True)
    output_file_path_1 = os.path.join(bond_length_folder_1, 'bond-length.txt')

    #Optimized molecular data path
    folder_path_2 = os.path.join(data_path_2, folder_name)
    bond_length_folder_2 = os.path.join(output_path_2, folder_name)
    os.makedirs(bond_length_folder_2, exist_ok=True)
    output_file_path_2 = os.path.join(bond_length_folder_2, 'bond-length-optimized.txt')

    for folder_path, output_file_path in [(folder_path_1, output_file_path_1), (folder_path_2, output_file_path_2)]:
        if not os.path.isdir(folder_path):
            print(f"Folder {folder_path} does not exist.")
            continue

        #Retrieve all. sdf files and sort them by file name
        sdf_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.sdf')], key=lambda x: int(x.split('.')[0]))

        with open(output_file_path, 'w') as output_file:
            for sdf_file in sdf_files:
                input_sdf_file_path = os.path.join(folder_path, sdf_file)
                try:
        
                    suppl = Chem.SDMolSupplier(input_sdf_file_path)
                except OSError as e:
                    if "Invalid input file" in str(e):
                        continue
                    else:
                        raise

                
                sdf_index = sdf_file.split('.')[0]
                for mol_index, mol in enumerate(suppl):
                    if mol is not None:
                        try:
                            bond_lengths = calculate_bond_lengths(mol)
                            output_line = "\t"
                            if sdf_index:
                                output_line = f"{sdf_index}\t"
                            for bond_key in all_bond_keys:
                                lengths = bond_lengths.get(bond_key, [])
                                lengths_str = ','.join(f"{length:.3f}" for length in lengths)
                                output_line += f"{bond_key}\t{lengths_str}\t"
                            output_file.write(output_line.strip() + '\n')
                        except Exception as e:
                            print(f"Error in {sdf_index} {mol_index}: {e}")
                            continue

#Merge files and extract column data
def calculate_Bond_length_merge(config):
    def merge_bond_length_files(base_path, output_filename, output_dir):
        folder_range = range(32)
        output_file_path = os.path.join(output_dir, output_filename)
        
        with open(output_file_path, 'w', encoding='utf-8') as outfile:
            for i in folder_range:
                folder_path = os.path.join(base_path, str(i))
                if os.path.exists(folder_path):
                    for file_name in os.listdir(folder_path):
                        if file_name.endswith('.txt'):
                            file_path = os.path.join(folder_path, file_name)
                            with open(file_path, 'r', encoding='utf-8') as infile:
                                outfile.write(infile.read() + '')

        if "optimized" in output_filename:
            target_suffix = "-optimized"
        else:
            target_suffix = ""
        extract_columns(output_file_path, output_dir, target_suffix)

    def extract_columns(file_path, output_dir, target_suffix=""):
        column_names = ['C-C', 'C=C', 'C≡C', 'C-N', 'C=N', 'C≡N', 'C-O', 'C=O', 'C:O']
        special_chars_map = {
            '≡': '#',
            '=': '=',
            ':': ':'
        }

        column_data = {name: [] for name in column_names}

        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                parts = line.strip().split('\t')
                for name in column_names:
                    try:
                        index = parts.index(name) + 1
                        if index < len(parts):
                            values = parts[index].split(',')
                        else:
                            values = []
                    except ValueError:
                        values = []
                    column_data[name].extend(values)

        os.makedirs(output_dir, exist_ok=True)
        for name in column_names:
            sanitized_name = name
            for char, replacement in special_chars_map.items():
                sanitized_name = sanitized_name.replace(char, replacement)
            
            output_file = os.path.join(output_dir, f'{sanitized_name}{target_suffix}.txt')
            with open(output_file, 'w', encoding='utf-8') as file:
                non_default_values = [value for value in column_data[name] if value]
                for value in non_default_values:
                    file.write(f"{value}\n")

    config_output_path = config['output']['output_path']
    base_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Bond-length')
    optimized_base_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Bond-length-optimized')
    output_dir = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Bond-length-merge')
    os.makedirs(output_dir, exist_ok=True)
    merge_bond_length_files(base_path, 'all-bond-length-merge.txt', output_dir)
    merge_bond_length_files(optimized_base_path, 'all-bond-length-optimized-merge.txt', output_dir)


def calculate_bond_length(config):
    try:
        config_output_path = config['output']['output_path']
        data_path_1 = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
        base_path = config_output_path
        data_path_2 = os.path.join(base_path, 'Structural-properties-metrics', '3D', 'Force-field-molecules')
        overall_output_path_1 = os.path.join(base_path, 'Structural-properties-metrics', '3D', 'Bond-length')
        overall_output_path_2 = os.path.join(base_path, 'Structural-properties-metrics', '3D', 'Bond-length-optimized')
        
        for folder_number in range(32):
            folder_name = str(folder_number)
            input_folder_1 = os.path.join(data_path_1, folder_name)
            input_folder_2 = os.path.join(data_path_2, folder_name)

            if not os.path.isdir(input_folder_1) and not os.path.isdir(input_folder_2):
                print(f"Neither input folder exists for {folder_name}. Skipping.")
                continue
            if os.path.isdir(input_folder_1):
                output_folder_1 = os.path.join(overall_output_path_1, folder_name)
                os.makedirs(output_folder_1, exist_ok=True)

            if os.path.isdir(input_folder_2):
                output_folder_2 = os.path.join(overall_output_path_2, folder_name)
                os.makedirs(output_folder_2, exist_ok=True)
            
            process(data_path_1, data_path_2, folder_name, overall_output_path_1, overall_output_path_2)

        calculate_Bond_length_merge(config)
    except Exception as e:
        print(f"Bond length calculation failed: {e}")

