import os
import numpy as np

def extract_values(file_path, keywords):
    values = {key: [] for key in keywords}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            parts = line.split()
            for i, part in enumerate(parts):
                for key in keywords:
                    if key in part:
                        try:
                            value = float(parts[i+1])
                            values[key].append(value)
                        except (ValueError, IndexError):
                            values[key].append(np.nan)  # Handle non-numeric or missing values gracefully
    return values

def extract_jsd_values(file_path):
    values = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            parts = line.split(':')
            if len(parts) > 1:
                try:
                    value = float(parts[1].strip())
                    values.append(value)
                except ValueError:
                    values.append(np.nan)  # Handle non-numeric values gracefully
    return values

def calculate_mean(values):
    if isinstance(values, dict):
        mean_values = {}
        for key, val_list in values.items():
            if len(val_list) > 0 and not all(np.isnan(val_list)):
                mean_values[key] = np.nanmean(val_list)
            else:
                mean_values[key] = 'NA'
        return mean_values
    elif isinstance(values, list):
        if len(values) > 0 and not all(np.isnan(values)):
            return np.nanmean(values)
        else:
            return 'NA'

def process_structural_metrics_2d(base_path, keywords):
    if not os.path.exists(base_path):
        #print(f"Directory {base_path} does not exist. Skipping 2D structure property metrics processing.")
        return None

    all_values = {key: [] for key in keywords}
    txt_files_found = False

    for file_name in os.listdir(base_path):
        if file_name.endswith('.txt'):
            txt_files_found = True
            file_path = os.path.join(base_path, file_name)
            file_values = extract_values(file_path, keywords)
            for key in keywords:
                all_values[key].extend(file_values[key])

    if not txt_files_found:
        return None

    return calculate_mean(all_values)

def process_structural_metrics_3d(base_paths):
    results = {}
    
    for key, path in base_paths.items():
        if not os.path.exists(path):
            #print(f"Directory {path} does not exist. Skipping 3D structure property metrics processing for {key}.")
            continue
        
        all_values = []
        txt_files_found = False

        for file_name in os.listdir(path):
            if file_name.endswith('.txt'):
                txt_files_found = True
                file_path = os.path.join(path, file_name)
                file_values = extract_jsd_values(file_path)
                all_values.extend(file_values)
        
        if not txt_files_found:
            continue
        
        results[key] = calculate_mean(all_values)

    return results

def calculate_summary(config):
    #Processing 2D structural attribute indicators
    output_path = config['output']['output_path']
    base_path_2d = os.path.join(output_path, "Structural-properties-metrics", "2D")
    output_file_2d = os.path.join(output_path, "Structural-properties-metrics", "2D-Structure-property-metrics.txt")
    keywords_2d = ['heavy-atom', 'chiral-center', 'Rings', 'Aromatic-Rings', 'Rotatable-Bonds', 'Fsp3']
    
    mean_values_2d = process_structural_metrics_2d(base_path_2d, keywords_2d)
    if mean_values_2d:
        with open(output_file_2d, 'w') as out_file:
            out_file.write('3D Structure property metrics\tMean Values\n')
            for key, mean_val in mean_values_2d.items():
                out_file.write(f"{key.replace('-', ' ')}\t{mean_val if mean_val == 'NA' else f'{mean_val:.6f}'}\n")

    #Processing 3D structural attribute indicators
    base_paths_3d = {
        'Bond_length_JSD': os.path.join(output_path, "Structural-properties-metrics", "3D", "JSD-of-bond-length"),
        'Bond_angles_JSD': os.path.join(output_path, "Structural-properties-metrics", "3D", "JSD-of-bond-angles"),
        'Dihedral_angles_JSD': os.path.join(output_path, "Structural-properties-metrics", "3D", "JSD-of-dihedral-angles"),
    }
    output_file_3d = os.path.join(output_path, "Structural-properties-metrics", "3D-Structure-property-metrics.txt")
    
    mean_values_3d = process_structural_metrics_3d(base_paths_3d)
    if mean_values_3d:
        with open(output_file_3d, 'w') as out_file:
            out_file.write('3D Structure property metrics\tMean Values\n')
            for key, mean_val in mean_values_3d.items():
                out_file.write(f"{key}\t{mean_val if mean_val == 'NA' else f'{mean_val:.6f}'}\n")
