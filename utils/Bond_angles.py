import os
import codecs
from rdkit.Chem import AllChem, SDMolSupplier
from rdkit import Chem
import re

def GetAngleDeg(conf, atom1, atom2, atom3):
    pt1 = conf.GetAtomPosition(atom1)
    pt2 = conf.GetAtomPosition(atom2)
    pt3 = conf.GetAtomPosition(atom3)
    v1 = [pt1.x - pt2.x, pt1.y - pt2.y, pt1.z - pt2.z]
    v2 = [pt3.x - pt2.x, pt3.y - pt2.y, pt3.z - pt2.z]
    dot_product = sum(a * b for a, b in zip(v1, v2))
    mag_v1 = sum(a**2 for a in v1)**0.5
    mag_v2 = sum(a**2 for a in v2)**0.5
    cos_theta = dot_product / (mag_v1 * mag_v2)
    angle_rad = abs(Chem.rdMolTransforms.GetAngleRad(conf, atom1, atom2, atom3))
    return angle_rad * 180 / 3.141592653589793

def get_mol_id(filename):
    match = re.search(r'([O\d+-]?\d+)\.sdf', filename)
    if match:
        return match.group(1)
    else:
        return None

def get_bond_angle_dict(mol, bond_angles_list):
    angle_dict = {}
    for bond_angles in bond_angles_list:
        substructure = Chem.MolFromSmiles(bond_angles)
        if substructure is not None:
            bond_pairs = mol.GetSubstructMatches(substructure)
            for pair in bond_pairs:
                try:
                    angle = GetAngleDeg(mol.GetConformer(), *pair)
                    assert mol.GetBondBetweenAtoms(pair[0], pair[1]) is not None
                    assert mol.GetBondBetweenAtoms(pair[2], pair[1]) is not None
                    if bond_angles not in angle_dict:
                        angle_dict[bond_angles] = []
                    angle_dict[bond_angles].append(angle)
                except Exception as e:
                    continue
    return angle_dict

def merge_bond_angles_files(base_path, output_filename, output_dir):
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
    column_names = ['CCC:', 'CSC:', 'CCO:', 'CNC:', 'OPO:', 'NCC:', 'CC=O:', 'COC:', 'CC=C:', 'OC=O:', 'NC=O:', 'CN=C:']
    
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

def calculate_bond_angles(config):
    try: 
        config_output_path = config['output']['output_path']
        data_path_1 = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
        data_path_2 = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Force-field-molecules')
        overall_output_path_1 = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Bond-angles')
        overall_output_path_2 = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Bond-angles-optimized')
        merge_output_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Bond-angles-merge')
        os.makedirs(overall_output_path_1, exist_ok=True)
        os.makedirs(overall_output_path_2, exist_ok=True)
        os.makedirs(merge_output_path, exist_ok=True)

        bond_angles_list = ['CCC', 'CSC', 'CCO', 'CNC', 'OPO', 'NCC', 'CC=O', 'COC', 'CC=C', 'OC=O', 'NC=O', 'CN=C']

        for data_path, overall_output_path in [(data_path_1, overall_output_path_1), (data_path_2, overall_output_path_2)]:
            for subdir in range(32): 
                subdir_path = os.path.join(data_path, str(subdir))
                if not os.path.isdir(subdir_path):
                    continue

                output_folder = os.path.join(overall_output_path, str(subdir))
                os.makedirs(output_folder, exist_ok=True)
                output_file_path = os.path.join(output_folder, 'bond-angles.txt')

                sdf_files = [f for f in os.listdir(subdir_path) if f.endswith(".sdf")]
                sorted_sdf_files = sorted(sdf_files, key=lambda x: int(get_mol_id(x)))

                with open(output_file_path, 'w', encoding='utf-8') as output_file:
                    for filename in sorted_sdf_files:
                        file_path = os.path.join(subdir_path, filename)
                        mol_id = get_mol_id(filename)
                        try:
                            suppl = SDMolSupplier(file_path)
                            for mol in suppl:
                                if mol is not None:
                                    try:
                                        mol = Chem.AddHs(mol)
                                        angle_dict = get_bond_angle_dict(mol, bond_angles_list)
                                        output_file.write(f"{mol_id}\t")
                                        for bond_angles in bond_angles_list:
                                            angles = angle_dict.get(bond_angles, ['NA'])
                                            output_file.write(f"{bond_angles}:\t{','.join(map(str, angles))}\t")
                                        output_file.write('\n')
                                    except Exception as e:
                                        continue
                        except OSError as e:
                            if "Invalid input file" in str(e):
                                continue  
                            else:
                                raise  

        merge_bond_angles_files(overall_output_path_1, 'all-bond-angles-merge.txt', merge_output_path)
        merge_bond_angles_files(overall_output_path_2, 'all-bond-angles-optimized-merge.txt', merge_output_path)

    except Exception as e:
        pass
