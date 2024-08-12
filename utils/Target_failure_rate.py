import os

def calculate_Target_failure_rate(config):
    
    input_path = config['data']['input_path']  
    output_base_path = config['output']['output_path']  
    output_dir = os.path.join(output_base_path, "Model-quality-metrics")  
    output_file = os.path.join(output_dir, "Target-failure-rate.txt")  

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = [] 
    empty_folder_count = 0  
    total_folders = 0  
    success = "success"
    fail = "fail"

    for i in range(32):
        folder = os.path.join(input_path, str(i))  
        
        if not os.path.exists(folder):
            continue  
        total_folders += 1

        sdf_files = [f for f in os.listdir(folder) if f.endswith('.sdf')]
        
        if not sdf_files:
            empty_folder_count += 1
            results.append(f"Target\t{i}, success/fail:\t{fail}\n")
        else:
            results.append(f"Target\t{i}, success/fail:\t{success}\n")
    if total_folders > 0:
        failure_rate = empty_folder_count / total_folders
    else:
        failure_rate = 0 

    results.append(f"Target failure rate mean:\t{failure_rate:.3f}\n")

    with open(output_file, 'w') as out_file:
        out_file.writelines(results)

