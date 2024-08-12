import os
from multiprocessing import Pool
from rdkit import Chem
import math

def process_folder(digit_folder_index, input_path, output_path):
    try:
        rmsd_results_main_folder = os.path.join(output_path, "Target-binding-metrics", "Glide-RDKit-RMSD-results")
        os.makedirs(rmsd_results_main_folder, exist_ok=True)
        rmsd_results_folder = os.path.join(rmsd_results_main_folder, str(digit_folder_index), "rmsd-results")
        os.makedirs(rmsd_results_folder, exist_ok=True)
        rmsd_folder = os.path.join(rmsd_results_main_folder, str(digit_folder_index), "rmsd-score")
        os.makedirs(rmsd_folder, exist_ok=True)
        min_rmsd_list = []
        input_folder = os.path.join(output_path, "Target-binding-metrics", "Glide-Docking-results", str(digit_folder_index))
        reference_folder = os.path.join(input_path, str(digit_folder_index))
        if not os.path.exists(input_folder) or not os.path.exists(reference_folder):
            
            return

        for input_file in os.listdir(input_folder):
            if input_file.endswith("_lib.sdf"):
                # Get the base file name of the input file (remove the suffix "_lib")
                input_file_base = os.path.splitext(input_file)[0].replace("_lib", "")
                input_file_path = os.path.join(input_folder, input_file)

                for reference_file in os.listdir(reference_folder):
                    if reference_file.endswith(".sdf"):
                        if input_file_base == os.path.splitext(reference_file)[0]:
                            output_file = os.path.join(rmsd_results_folder, f"rmsd-results-{input_file_base}-{reference_file}.txt")
                            reffile = os.path.join(reference_folder, reference_file)
                            probfile = input_file_path

                            refmol_supplier = Chem.SDMolSupplier(reffile)
                            probmol_supplier = Chem.SDMolSupplier(probfile)

                            if refmol_supplier is None or probmol_supplier is None:
                                #print(f"Error reading {reffile} or {probfile}")
                                continue

                            refmol = refmol_supplier[0]
                            probmol = probmol_supplier[0]

                            if refmol is None or probmol is None:
                                #print(f"Error: Molecule could not be read from {reffile} or {probfile}")
                                continue

                            refmol = Chem.RemoveHs(refmol)
                            probmol = Chem.RemoveHs(probmol)

                            refatomnum = refmol.GetNumAtoms()
                            probatomnum = probmol.GetNumAtoms()

                            if refatomnum != probatomnum:
                                #print(f"Atom number mismatch between {reffile} and {probfile}")
                                continue

                            refconf = refmol.GetConformer()
                            probconf = probmol.GetConformer()

                            d2_list = []
                            for i in range(refatomnum):
                                ref_atom_coord = refconf.GetAtomPosition(i)
                                prob_atom_coord = probconf.GetAtomPosition(i)
                                d_x = ref_atom_coord.x - prob_atom_coord.x
                                d_y = ref_atom_coord.y - prob_atom_coord.y
                                d_z = ref_atom_coord.z - prob_atom_coord.z
                                d = math.sqrt(d_x ** 2 + d_y ** 2 + d_z ** 2)
                                d2_list.append(d ** 2)
                            rmsd_value = math.sqrt(sum(d2_list) / refatomnum)
                            min_rmsd_list.append(f"RMSD for\t{input_file}\tand\t{reference_file}:\t{rmsd_value:.3f}")

                            with open(output_file, "w") as file:
                                file.write(f"RMSD:\t{rmsd_value:.3f}\n")
                            #print(f"Wrote RMSD to {output_file}")

        # Summarize all RMSD results and write them into a summary file
        rmsd_summary_file = os.path.join(rmsd_folder, "rmsd-summary.txt")
        with open(rmsd_summary_file, "w") as summary_file:
            for rmsd_entry in min_rmsd_list:
                summary_file.write(rmsd_entry + '\n')
            # print(f"RMSD value {rmsd_summary_file} ")
    except Exception as e:
        pass

def calculate_rdkit_rmsd(config):
    try:
        output_path = config['output']['output_path']
        input_path = os.path.join(output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
        digit_folders = list(range(32))

        # Filter out digit_folders that actually exist in both input and output paths
        existing_folders = []
        for digit_folder in digit_folders:
            input_folder = os.path.join(output_path, "Target-binding-metrics", "Glide-Docking-results", str(digit_folder))
            reference_folder = os.path.join(input_path, str(digit_folder))
            if os.path.exists(input_folder) and os.path.exists(reference_folder):
                existing_folders.append(digit_folder)

        with Pool() as pool:
            pool.starmap(process_folder, [(digit_folder, input_path, output_path) for digit_folder in existing_folders])
    except Exception as e:
        pass


