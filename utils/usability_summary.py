import os

# Function to count lines in a file
def count_lines(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            return sum(1 for _ in file)  # Count the number of lines in the file
    return 0

# Function to calculate and summarize the usability ratio
def calculate_usability_summary(config):
    """
    Calculate and summarize the usability ratio.
    """
    config_output_path = config['output']['output_path']
    input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics')
    Uniqueness_base_path = os.path.join(input_path, 'Uniqueness')
    Usability_base_path = os.path.join(input_path, 'Usability')
    summary_file_path = os.path.join(Usability_base_path, 'Usability_summary.txt')

    summary_data = []
    total_ratio = 0
    valid_folder_count = 0

    for folder_number in range(32):  # Iterate through folders 0 to 31
        uniqueness_folder = os.path.join(Uniqueness_base_path, str(folder_number))
        usability_folder = os.path.join(Usability_base_path, str(folder_number))

        if not os.path.exists(uniqueness_folder) or not os.path.exists(usability_folder):
            continue  # Skip this folder if it does not exist

        unique_file_path = os.path.join(uniqueness_folder, 'unique.txt')
        usability_file_path = os.path.join(usability_folder, 'usability.txt')

        unique_lines = count_lines(unique_file_path)
        usability_lines = count_lines(usability_file_path)

        ratio = unique_lines / usability_lines if usability_lines > 0 else 0  # Calculate usability ratio
        summary_data.append((folder_number, unique_lines, usability_lines, ratio))
        total_ratio += ratio
        valid_folder_count += 1

    average_ratio = total_ratio / valid_folder_count if valid_folder_count > 0 else 0  # Calculate average usability ratio

    with open(summary_file_path, 'w') as summary_file:  # Write to summary file
        summary_file.write("Targets\tUnique_Count\tUsability_Count\tUsability\n")
        for folder_number, unique_lines, usability_lines, ratio in summary_data:
            summary_file.write(f"{folder_number}\t{unique_lines}\t{usability_lines}\t{ratio:.2f}\n")
        summary_file.write(f"Mean of Usability:\t{average_ratio:.3f}\n")  # Write average usability ratio
