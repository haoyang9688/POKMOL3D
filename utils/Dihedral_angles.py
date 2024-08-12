import os
import codecs
from rdkit import Chem
from rdkit.Chem import AllChem

def get_bond_symbol(bond):
    """
    Return the symbol representation of a bond
    """
    a0 = bond.GetBeginAtom().GetSymbol()
    a1 = bond.GetEndAtom().GetSymbol()
    b = str(int(bond.GetBondType()))  # single: 1, double: 2, triple: 3, aromatic: 12
    return ''.join([a0, b, a1])

def get_triple_bonds(mol):
    """
    Get all the bond triplets in a molecule
    """
    valid_triple_bonds = []
    for idx_bond, bond in enumerate(mol.GetBonds()):
        idx_begin_atom = bond.GetBeginAtomIdx()
        idx_end_atom = bond.GetEndAtomIdx()
        begin_atom = mol.GetAtomWithIdx(idx_begin_atom)
        end_atom = mol.GetAtomWithIdx(idx_end_atom)
        begin_bonds = begin_atom.GetBonds()
        valid_left_bonds = []
        for begin_bond in begin_bonds:
            if begin_bond.GetIdx() == idx_bond:
                continue
            else:
                valid_left_bonds.append(begin_bond)
        if len(valid_left_bonds) == 0:
            continue

        end_bonds = end_atom.GetBonds()
        for end_bond in end_bonds:
            if end_bond.GetIdx() == idx_bond:
                continue
            else:
                for left_bond in valid_left_bonds:
                    valid_triple_bonds.append([left_bond, bond, end_bond])
    return valid_triple_bonds

def Find_bond_triplets(mol, bonds_ref_sym_list):
    """
    Find bond triplets (defined by bonds_ref_sym_list) in mol and return the dihedral angles of each triplet.
    """
    bonds_list = get_triple_bonds(mol)
    angles_dict = {sym: [] for sym in bonds_ref_sym_list}

    for bonds in bonds_list:
        sym = '-'.join([get_bond_symbol(b) for b in bonds])
        sym1 = '-'.join([get_bond_symbol(b) for b in bonds][::-1])

        if sym in angles_dict or sym1 in angles_dict:
            if sym1 in angles_dict:
                bonds = bonds[::-1]
                sym = sym1

            bond0 = bonds[0]
            atom0 = bond0.GetBeginAtomIdx()
            atom1 = bond0.GetEndAtomIdx()

            bond1 = bonds[1]
            atom1_0 = bond1.GetBeginAtomIdx()
            atom1_1 = bond1.GetEndAtomIdx()
            if atom0 == atom1_0:
                i, j, k = atom1, atom0, atom1_1
            elif atom0 == atom1_1:
                i, j, k = atom1, atom0, atom1_0
            elif atom1 == atom1_0:
                i, j, k = atom0, atom1, atom1_1
            elif atom1 == atom1_1:
                i, j, k = atom0, atom1, atom1_0

            bond2 = bonds[2]
            atom2_0 = bond2.GetBeginAtomIdx()
            atom2_1 = bond2.GetEndAtomIdx()
            if atom2_0 == k:
                l = atom2_1
            elif atom2_1 == k:
                l = atom2_0

            angle = Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(), i, j, k, l)
            angles_dict[sym].append(angle)

    return angles_dict

def calculate_dihedral_angles(config):
    config_output_path = config['output']['output_path']
    data_path_1 = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
    data_path_2 = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Force-field-molecules')
    overall_output_path_1 = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Dihedral-angles')
    overall_output_path_2 = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Dihedral-angles-optimized')
    merge_output_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Dihedral-angles-merge')
    os.makedirs(overall_output_path_1, exist_ok=True)
    os.makedirs(overall_output_path_2, exist_ok=True)
    os.makedirs(merge_output_path, exist_ok=True)

    bond_symbols_list = ['C1C-C1C-C1C', 'C12C-C12C-C12C', 'C1C-C1C-C1O', 'O1C-C1C-C1O', 'C1C-C12C-C12C', 'C1C-C2C-C1C']

    for data_path, overall_output_path in [(data_path_1, overall_output_path_1), (data_path_2, overall_output_path_2)]:
        for foldername in os.listdir(data_path):
            subfolder = os.path.join(data_path, foldername)
            if os.path.isdir(subfolder) and foldername.isdigit() and int(foldername) in range(32):
                output_folder = os.path.join(overall_output_path, foldername)
                os.makedirs(output_folder, exist_ok=True)
                output_file_path = os.path.join(output_folder, 'dihedral-angles.txt')

                sdf_files = [f for f in os.listdir(subfolder) if f.endswith(".sdf")]
                sdf_files.sort(key=lambda x: int(os.path.splitext(x)[0]))  # Sort by number
                
                with open(output_file_path, 'w', encoding='utf-8') as output_file:
                    for filename in sdf_files:
                        molecule_number = os.path.splitext(filename)[0]  # Extract file names without extensions
                        file_path = os.path.join(subfolder, filename)
                        
                        try:
                            suppl = Chem.SDMolSupplier(file_path)
                            if not suppl:
                                continue
                        except OSError as e:
                            #print(f"Error: {e}")
                            continue

                        for mol in suppl:
                            if mol is None:
                                continue
                            mol = Chem.AddHs(mol)
                            
                            angles_dict = Find_bond_triplets(mol, bond_symbols_list)
                            mol = Chem.RemoveHs(mol)
                            output_file.write(f"{molecule_number}\t")
                            for bond_sym in bond_symbols_list:
                                if bond_sym in angles_dict:
                                    angles = angles_dict[bond_sym]
                                    output_file.write(f"{bond_sym}\t{angles[0]:.3f}\t" if angles else f"{bond_sym}\tNA\t")
                                else:
                                    output_file.write(f"{bond_sym}\tNA\t")
                            output_file.write("\n")

    merge_dihedral_angles_files(overall_output_path_1, 'all-dihedral-angles-merge.txt', merge_output_path)
    merge_dihedral_angles_files(overall_output_path_2, 'all-dihedral-angles-optimized-merge.txt', merge_output_path)

def merge_dihedral_angles_files(base_path, output_filename, output_dir):
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
                            outfile.write(infile.read() + '')  # Do not add line breaks

    if "optimized" in output_filename:
        target_suffix = "-optimized"
    else:
        target_suffix = ""
    extract_columns(output_file_path, output_dir, target_suffix)

def extract_columns(file_path, output_dir, target_suffix=""):
    column_names = ['C1C-C1C-C1C', 'C12C-C12C-C12C', 'C1C-C1C-C1O', 'O1C-C1C-C1O', 'C1C-C12C-C12C', 'C1C-C2C-C1C']
    
    column_data = {name: [] for name in column_names}

    with codecs.open(file_path, 'r', encoding='utf-8') as file:
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
                cleaned_values = [value.strip() for value in values if value.strip() != 'NA']
                column_data[name].extend(cleaned_values)

    os.makedirs(output_dir, exist_ok=True)
    for name in column_names:
        output_file_name = f"{name.strip(':')}{target_suffix}.txt"
        output_file_path = os.path.join(output_dir, output_file_name)
        with codecs.open(output_file_path, 'w', encoding='utf-8') as file:
            for value in column_data[name]:
                file.write(f"{value}\n")
