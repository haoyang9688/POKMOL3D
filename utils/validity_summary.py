import os

# Function to count lines in a file
def count_lines(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            return sum(1 for _ in file)  # Count the number of lines in the file
    return 0

# Function to count files with a specific extension in a folder
def count_files(folder_path, extension):
    if os.path.exists(folder_path):
        return len([f for f in os.listdir(folder_path) if f.endswith(extension)])  # Count the number of files with the given extension
    return 0

# Function to calculate validity summary
def calculate_validity_summary(config):
    base_path = config['output']['output_path']
    validity_path = os.path.join(base_path, 'General-molecular-quality-metrics', 'Validity')
    sdf_path = config['data']['input_path']
    summary_file_path = os.path.join(validity_path, 'Validity-summary.txt')

    summary_data = []
    total_ratio = 0
    valid_folders_count = 0  # Counter for the actual number of valid folders

    for folder_number in range(32):  # Iterate through folders 0 to 31
        validity_folder = os.path.join(validity_path, str(folder_number))
        sdf_folder = os.path.join(sdf_path, str(folder_number))
        
        # Check if both the validity folder and sdf folder exist, if not, skip to the next folder
        if not os.path.exists(validity_folder) or not os.path.exists(sdf_folder):
            continue

        validity_file_path = os.path.join(validity_folder, 'validity.txt')
        sdf_count = count_files(sdf_folder, '.sdf')
        validity_lines = count_lines(validity_file_path)

        ratio = validity_lines / sdf_count if sdf_count > 0 else 0  # Calculate validity ratio
        summary_data.append((folder_number, sdf_count, validity_lines, ratio))
        total_ratio += ratio
        valid_folders_count += 1  # Increment the counter for valid folders

    average_ratio = total_ratio / valid_folders_count if valid_folders_count > 0 else 0  # Calculate average validity ratio

    with open(summary_file_path, 'w') as summary_file:  # Write to summary file
        summary_file.write("Targets\tSDF_Count\tValidity_count\tValidity\n")
        for folder_number, sdf_count, validity_lines, ratio in summary_data:
            summary_file.write(f"{folder_number}\t{sdf_count}\t{validity_lines}\t{ratio:.2f}\n")
        summary_file.write(f"Mean of validity:\t{average_ratio:.3f}\n")  # Write average validity ratio


