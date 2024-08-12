import os
import subprocess
from multiprocessing import Pool
def process_folder(digit_folder_index, input_path, output_path, rmsd_script_path):
    rmsd_results_main_folder = os.path.join(output_path, "Target-binding-metrics", "Vina-Openeye-RMSD-results")
    os.makedirs(rmsd_results_main_folder, exist_ok=True)

    rmsd_results_folder = os.path.join(rmsd_results_main_folder, str(digit_folder_index), "rmsd-results")
    os.makedirs(rmsd_results_folder, exist_ok=True)
    
    rmsd_folder = os.path.join(rmsd_results_main_folder, str(digit_folder_index), "rmsd-score")
    os.makedirs(rmsd_folder, exist_ok=True)

    rmsd_list = []
    input_folder = os.path.join(output_path, "Target-binding-metrics", "Vina-Docking-results", str(digit_folder_index), "docking-poses")
    reference_folder = os.path.join(input_path, str(digit_folder_index))

    #Debugging information: Print input folder and reference folder paths
    #print(f"Processing folder index: {digit_folder_index}")
    #print(f"Input folder: {input_folder}")
    #print(f"Reference folder: {reference_folder}")

    if not os.path.exists(input_folder):
        #print(f"Input folder does not exist: {input_folder}")
        return

    if not os.path.exists(reference_folder):
        #print(f"Reference folder does not exist: {reference_folder}")
        return

    #Traverse the files in the input folder
    for input_file in os.listdir(input_folder):
        if input_file.endswith(".sdf"):
            input_file_base = os.path.splitext(input_file)[0]
            input_file_path = os.path.join(input_folder, input_file)
            reference_file = f"{input_file_base}.sdf"
            reference_file_path = os.path.join(reference_folder, reference_file)

            if os.path.exists(reference_file_path):
                output_file = os.path.join(rmsd_results_folder, f"rmsd-results-{input_file_base}.sdf-{reference_file}.txt")
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                
                command = [
                    "python", rmsd_script_path,
                    "-ref", reference_file_path,
                    "-in", input_file_path,
                    "-out", output_file
                ]

                try:
                    #Run RMSD calculation script
                    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    if result.returncode == 0:
                        print(f"RMSD calculation succeeded for {input_file} and {reference_file}")
                    else:
                        print(f"RMSD calculation failed for {input_file} and {reference_file}")
                        print(result.stderr)
                except subprocess.CalledProcessError as e:
                    print(f"Error during RMSD calculation for {input_file} and {reference_file}: {e}")
                    continue

                #Read the result file and extract the RMSD value
                try:
                    with open(output_file, "r") as results_file:
                        lines = results_file.readlines()
                except FileNotFoundError:
                    print(f"Output file not found: {output_file}")
                    continue

                rmsd = None
                for line in lines:
                    if "RMSD-0:" in line:
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            rmsd_str = parts[-1].strip()
                            try:
                                rmsd_value = float(rmsd_str)
                                if rmsd is None or rmsd_value < rmsd:
                                    rmsd = rmsd_value
                            except ValueError:
                                pass

                #Add the extracted RMSD values to the list
                if rmsd is not None:
                    rmsd_list.append(f"RMSD for\t{input_file_base}.sdf_{reference_file}:\t{rmsd}")
            else:
                print(f"Reference file does not exist: {reference_file_path}")

    rmsd_summary_file = os.path.join(rmsd_folder, "rmsd-summary.txt")
    os.makedirs(os.path.dirname(rmsd_summary_file), exist_ok=True)
    
    with open(rmsd_summary_file, "w") as summary_file:
        for rmsd_entry in rmsd_list:
            summary_file.write(rmsd_entry + '\n')

def calculate_rmsd(config):
    
    output_path = config['output']['output_path']
    input_path = os.path.join(output_path, 'General-molecular-quality-metrics', 'Uniqueness-sdf')
    POKMOL3D_path = config['POKMOL3D_path']
    rmsd_script_path = os.path.join(POKMOL3D_path, "utils", "rmsd.py")
    digit_folders = list(range(32))
    # Filter out non-existent folders before processing
    existing_folders = [
        digit_folder for digit_folder in digit_folders
        if os.path.exists(os.path.join(output_path, "Target-binding-metrics", "Vina-Docking-results", str(digit_folder))) and
           os.path.exists(os.path.join(input_path, str(digit_folder)))
    ]

    with Pool() as pool:
        pool.starmap(process_folder, [(digit_folder, input_path, output_path, rmsd_script_path) for digit_folder in existing_folders])

