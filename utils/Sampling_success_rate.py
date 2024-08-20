import os

def calculate_Sampling_success_rate(config):
    input_path = config['data']['input_path']
    output_base_path = config['output']['output_path']
    output_dir = os.path.join(output_base_path, "Model-quality-metrics")
    output_file = os.path.join(output_dir, "Sampling-success-rate.txt")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = []
    success_folder_count = 0
    total_folders = 0 

    for i in range(32):
        folder = os.path.join(input_path, str(i))
        
        if not os.path.exists(folder):
            continue
        
        total_folders += 1  
    
        sdf_files = [f for f in os.listdir(folder) if f.endswith('.sdf')]
        success_rate = 1 if len(sdf_files) >= 2000 else 0
        success_folder_count += success_rate
        
        results.append(f"Target\t{i}\tSampling_success_rate:\t{success_rate}\n")

    if total_folders > 0:
        average_success_rate = success_folder_count / total_folders
    else:
        average_success_rate = 0  

    results.append(f"Sampling success rate mean:\t{average_success_rate:.3f}\n")

    with open(output_file, 'w') as out_file:
        out_file.writelines(results)




