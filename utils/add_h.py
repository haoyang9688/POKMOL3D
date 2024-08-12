import os
from pymol import cmd
import multiprocessing

def hydrogenate_with_pymol(input_sdf, output_folder):
    """
    Use PyMOL to hydrogenate the input SDF file and save the result to the specified output folder。
    
    input_sdf (str): The input SDF file path。
    output_folder (str): Output folder path。
    """
    #Load SDF file into PyMOL
    cmd.load(input_sdf, 'mol')
    
    #Add hydrogen atoms to all atoms
    cmd.h_add(selection='all')
    
    #Save the processed molecules in SDF format to the output folder
    output_sdf = os.path.join(output_folder, os.path.basename(input_sdf))
    cmd.save(output_sdf, 'mol', format='sdf')
    
    #Delete loaded molecular objects
    cmd.delete('mol')

def process_sdf_file(args):
    """
    A function for processing a single SDF file, used for multiprocessing parallel processing。
    parameter:
    ARGS (tuple): A tuple containing the input SDF file path and output folder path。
    """
    input_sdf, output_folder = args
    hydrogenate_with_pymol(input_sdf, output_folder)

def calculate_add_h_main(config):
   
    output_root_folder = config['output']['output_path']
    input_root_folder = os.path.join(output_root_folder, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
    number_folders = [f for f in os.listdir(input_root_folder) if os.path.isdir(os.path.join(input_root_folder, f))]

    tasks = []

    for number_folder in number_folders:
        number_folder_path = os.path.join(input_root_folder, number_folder)
        output_folder = os.path.join(output_root_folder, "Target-binding-metrics", "Prepared-add-h-ligands", number_folder)
        os.makedirs(output_folder, exist_ok=True)
        sdf_files = [f for f in os.listdir(number_folder_path) if f.endswith('.sdf')]
        
        for sdf_file in sdf_files:
            input_sdf = os.path.join(number_folder_path, sdf_file)
            tasks.append((input_sdf, output_folder))

    with multiprocessing.Pool() as pool:
        pool.map(process_sdf_file, tasks)
    
    #print("Hydrogenation completed.")

