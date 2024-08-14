import os
from rdkit import Chem
from rdkit.Chem import AllChem

def calculate_MMFF_preparation(config):
    config_input_path = config['data']['input_path']
    config_output_path = config['output']['output_path']
    sdf_output_base_dir = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Force-field-molecules')
    os.makedirs(sdf_output_base_dir, exist_ok=True)
    
    total_molecules = 0
    embedded_molecules = 0
    optimized_molecules = 0

    
    for i in range(32): 
        folder_name = os.path.join(config_input_path, str(i))  
        if os.path.exists(folder_name): 
            sdf_output_dir = os.path.join(sdf_output_base_dir, str(i))
            
            
            if not os.path.exists(sdf_output_dir):
                os.makedirs(sdf_output_dir)
            
            for file_name in os.listdir(folder_name):  
                if file_name.endswith('.sdf'): 
                    sdf_input = os.path.join(folder_name, file_name)  
                    sdf_output = os.path.join(sdf_output_dir, f'{file_name}')
                    
                    #print(f": {sdf_input}")
                    #print(f": {sdf_output}")
                    supplier = Chem.SDMolSupplier(sdf_input)
                  
                    if not supplier:
                        #print(f"error: {sdf_input}")
                        continue

                    writer = Chem.SDWriter(sdf_output)

                    for mol in supplier:
                        if mol is None:
                            #print("Skipping invalid molecules")
                            continue

                        total_molecules += 1

                        mol = Chem.AddHs(mol)
                        try:
                            embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                            if embed_result == -1:
                                #print("Unable to generate 3D conformation")
                                continue
                        except Exception as e:
                            #print(f"error: {e}")
                            continue

                        embedded_molecules += 1

                        try:
                            optimize_result = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=200)
                            if optimize_result != 0:
                                #print("error")
                                continue
                        except Exception as e:
                            #print(f"error: {e}")
                            continue

                        optimized_molecules += 1
                        writer.write(mol)

                    writer.close()
                    #print(f"The optimized molecule has been written into: {sdf_output}")
        else:
            #print(f"skip：{folder_name}")
            pass


