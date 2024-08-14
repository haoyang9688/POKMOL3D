import os

# Function to count lines in a file
def count_lines(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            return sum(1 for _ in file)  # Count the number of lines in the file
    return 0

# Function to calculate and summarize the uniqueness and effectiveness ratio
def calculate_uniqueness_summary(config):
    config_output_path = config['output']['output_path']
    input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics')
    Uniqueness_base_path = os.path.join(input_path, 'Uniqueness')
    Validity_base_path = os.path.join(input_path, 'Validity')
    summary_file_path = os.path.join(Uniqueness_base_path, 'Uniqueness-summary.txt')

    summary_data = []
    total_ratio = 0
    actual_folder_count = 0

    for folder_number in range(32):  # Iterate through folders 0 to 31
        uniqueness_folder = os.path.join(Uniqueness_base_path, str(folder_number))
        validity_folder = os.path.join(Validity_base_path, str(folder_number))

        unique_file_path = os.path.join(uniqueness_folder, 'unique.txt')
        validity_file_path = os.path.join(validity_folder, 'validity.txt')

        if os.path.exists(uniqueness_folder) and os.path.exists(validity_folder):
            unique_lines = count_lines(unique_file_path)
            validity_lines = count_lines(validity_file_path)

            if unique_lines > 0 or validity_lines > 0:  # Ensure that there is data to calculate
                ratio = unique_lines / validity_lines if validity_lines > 0 else 0  # Calculate uniqueness ratio
                summary_data.append((folder_number, unique_lines, validity_lines, ratio))
                total_ratio += ratio
                actual_folder_count += 1

    average_ratio = total_ratio / actual_folder_count if actual_folder_count > 0 else 0  # Calculate average uniqueness ratio

    if summary_data:  # Only write if there is data to write
        with open(summary_file_path, 'w') as summary_file:  # Write to summary file
            summary_file.write("Targets\tUnique_Count\tValidity_count\tUniqueness\n")
            for folder_number, unique_lines, validity_lines, ratio in summary_data:
                summary_file.write(f"{folder_number}\t{unique_lines}\t{validity_lines}\t{ratio:.3f}\n")
            summary_file.write(f"Mean of Uniqueness:\t{average_ratio:.3f}\n")  # Write average uniqueness ratio


