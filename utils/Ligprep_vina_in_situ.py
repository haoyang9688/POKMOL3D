import os
import subprocess
from multiprocessing import Pool
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from rdkit import Chem

def process_sdf_file(input_sdf_file, output_pdbqt_file):
    for mol in Chem.SDMolSupplier(input_sdf_file, removeHs=False):
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(output_pdbqt_file, "w") as file:
                    file.write(pdbqt_string)

def process_folder(folder_path, output_base_folder):
    sdf_files = [f for f in os.listdir(folder_path) if f.endswith(".sdf") and os.path.isfile(os.path.join(folder_path, f))]
    docking_output_pdbqt_folder = os.path.join(output_base_folder, os.path.basename(folder_path), "Ligprep-output-pdbqt")
    os.makedirs(docking_output_pdbqt_folder, exist_ok=True)
    for sdf_file in sdf_files:
        input_sdf_file = os.path.join(folder_path, sdf_file)
        output_pdbqt_file = os.path.join(docking_output_pdbqt_folder, os.path.splitext(sdf_file)[0] + ".pdbqt")
        process_sdf_file(input_sdf_file, output_pdbqt_file)

def vina_insitu_ligprep(config):
    config_output_path = config['output']['output_path']
    main_directory = os.path.join(config_output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
    config_redocking_output_path = os.path.join(config_output_path, 'Target-binding-metrics', 'Prepared-add-h-ligands')
    os.makedirs(config_redocking_output_path, exist_ok=True)

    folders = [os.path.join(main_directory, str(i)) for i in range(32) if os.path.isdir(os.path.join(main_directory, str(i)))]

    with Pool(processes=100) as pool:
        pool.starmap(process_folder, [(folder, config_redocking_output_path) for folder in folders])

