import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

OUTPUT_FILE_NAME = '2D_properties.txt'

def structural_properties(file_name, smiles, output_folder):
    """
    Calculate the topological structure characteristics of molecules and save the results to an output file.

    :param file_name: File name of the SDF file
    :param smiles: SMILES
    :param output_folder: Output file path
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Failed to load molecule from {file_name}")
            return

        # Calculate various topological structural features
        num_atoms = mol.GetNumAtoms()
        num_non_hydrogen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        num_chiral_carbons = sum(1 for idx, _ in chiral_centers if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        fraction_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)

        # Format the properties string
        properties_str = (f"{file_name}\theavy-atom:\t{num_non_hydrogen_atoms}\tchiral-center:\t{num_chiral_carbons}\t"
                          f"Rings:\t{num_rings}\tAromatic-Rings:\t{num_aromatic_rings}\t"
                          f"Rotatable-Bonds:\t{num_rotatable_bonds}\tFsp3:\t{fraction_csp3:.3f}")
        
        # Write the properties string to the output file
        output_path = os.path.join(output_folder, OUTPUT_FILE_NAME)
        with open(output_path, 'a') as result_file:
            result_file.write(properties_str + '\n')
    except Exception as e:
        print(f"Error processing molecule from {file_name}: {e}")

def process(data_path, folder_name, overall_output_path):
    """
    Process all .sdf files in the specified folder.

    :param data_path: Root folder path
    :param folder_name: Name of the subfolder to be processed
    :param overall_output_path: Output file path
    """
    folder_path = os.path.join(data_path, folder_name)

    if not os.path.isdir(folder_path):
        print(f"Folder {folder_path} does not exist, skipping.")
        return

    output_folder = os.path.join(overall_output_path, folder_name)
    os.makedirs(output_folder, exist_ok=True)

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".sdf"):
            file_path = os.path.join(folder_path, file_name)
            try:
                suppl = Chem.SDMolSupplier(file_path)
                for idx, mol in enumerate(suppl):
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        structural_properties(file_name, smiles, output_folder)
                    else:
                        print(f"Failed to parse molecule {idx} in file {file_name}")
            except Exception as e:
                print(f"Error processing file {file_name}: {e}")

def calculate_structural_properties_2D(config):
    """
    Calculate topological structure characteristics for molecules in multiple subfolders.

    :param config: Configuration parameters, including input and output paths
    """
    try:
        config_output_path = config['output']['output_path']
        input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
        overall_output_path = os.path.join(config_output_path, 'Structural-properties-metrics', '2D')
        os.makedirs(overall_output_path, exist_ok=True)

        for folder_number in range(32):
            folder_name = str(folder_number)
            input_folder_path = os.path.join(input_path, folder_name)
            
            if os.path.isdir(input_folder_path):  # Check if the folder exists
                output_folder = os.path.join(overall_output_path, folder_name)
                os.makedirs(output_folder, exist_ok=True)  # Create output folder only if the input folder exists
                process(input_path, folder_name, overall_output_path)
            else:
                print(f"Input folder {input_folder_path} does not exist, skipping.")
    except Exception as e:
        print(f"structural-properties-2D failed: {e}")
