import numpy as np
from scipy.spatial.distance import jensenshannon
import os
from scipy.stats import gaussian_kde

def calculate_js_divergence(data1, data2, num_points=100, epsilon=1e-6):
    kde_P = gaussian_kde(data1)
    kde_Q = gaussian_kde(data2)
    x_eval = np.linspace(min(data1.min(), data2.min()), max(data1.max(), data2.max()), num=num_points)
    P = kde_P(x_eval) + epsilon
    Q = kde_Q(x_eval) + epsilon
    P /= np.sum(P)
    Q /= np.sum(Q)
    js_divergence = jensenshannon(P, Q)
    return js_divergence

def read_data(file_path):
    data = np.genfromtxt(file_path)
    
    data = data[~np.isnan(data)]
    return data


def format_output(file_pair, js_divergence):
    file_name = os.path.splitext(file_pair[1])[0]
    output_line = f"{file_name}\tJS Divergence:\t{js_divergence:.3f}"
    return output_line + "\n"


def calculate_dajs_divergence(config):
    config_output_path = config['output']['output_path']
    folder_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Dihedral-angles-merge')

    
    file_pairs = [
        ('C1C-C1C-C1C-optimized.txt', 'C1C-C1C-C1C.txt'),
        ('C12C-C12C-C12C-optimized.txt', 'C12C-C12C-C12C.txt'),
        ('C1C-C1C-C1O-optimized.txt', 'C1C-C1C-C1O.txt'),
        ('O1C-C1C-C1O-optimized.txt', 'O1C-C1C-C1O.txt'),
        ('C1C-C12C-C12C-optimized.txt', 'C1C-C12C-C12C.txt'),
        ('C1C-C2C-C1C-optimized.txt', 'C1C-C2C-C1C.txt')
    ]

    output_lines = []
    for file_pair in file_pairs:
        data1 = read_data(os.path.join(folder_path, file_pair[0]))
        data2 = read_data(os.path.join(folder_path, file_pair[1]))
        js_divergence = calculate_js_divergence(data1, data2)  
        output_lines.append(format_output(file_pair, js_divergence))

    output_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'JSD-of-dihedral-angles')
    os.makedirs(output_path, exist_ok=True)

    output_file_path = os.path.join(output_path, 'JS-divergence-result.txt')
    with open(output_file_path, 'w') as output_file:
        for line in output_lines:
            output_file.writelines(line)
