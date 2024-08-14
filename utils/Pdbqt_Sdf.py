import os
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
from rdkit import Chem

def convert_pdbqt_to_sdf(input_file, output_file):
    try:
        pdbqt_mol = PDBQTMolecule.from_file(input_file, skip_typing=True)
        rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        writer = Chem.SDWriter(output_file)
        for rdkitmol in rdkitmol_list:
            writer.write(rdkitmol)
        writer.close()
        #print(f"Successfully converted {input_file} to {output_file}")
    except Exception as e:
        pass

def convert_pdbqt_sdf(config):
    base_dir = os.path.join(config['output']['output_path'], 'Target-binding-metrics', 'Vina-Docking-results')
    
    for i in range(32):
        folder_path = os.path.join(base_dir, str(i), 'docking-poses')

        if os.path.exists(folder_path):
            for file_name in os.listdir(folder_path):
                if file_name.endswith('.pdbqt'):
                    input_file = os.path.join(folder_path, file_name)
                    if file_name.startswith('output-'):
                        output_file_name = file_name.replace('output-', '').replace('.pdbqt', '.sdf')
                    else:
                        output_file_name = file_name.replace('.pdbqt', '.sdf')

                    output_file = os.path.join(folder_path, output_file_name)

                    convert_pdbqt_to_sdf(input_file, output_file)
        else:
            print(f"pass：{folder_path}")
