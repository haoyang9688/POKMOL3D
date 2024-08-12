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


def calculate_bajs_divergence(config):
    config_output_path = config['output']['output_path']
    folder_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'Bond-angles-merge')

    file_pairs = [
        ('CCC-optimized.txt', 'CCC.txt'),
        ('CSC-optimized.txt', 'CSC.txt'),
        ('CCO-optimized.txt', 'CCO.txt'),
        ('CNC-optimized.txt', 'CNC.txt'),
        ('NCC-optimized.txt', 'NCC.txt'),
        ('CC=O-optimized.txt', 'CC=O.txt'),
        ('COC-optimized.txt', 'COC.txt'),
        ('CC=C-optimized.txt', 'CC=C.txt'),
        ('OC=O-optimized.txt', 'OC=O.txt'),
        ('NC=O-optimized.txt', 'NC=O.txt'),
        ('CN=C-optimized.txt', 'CN=C.txt')
    ]

    output_lines = []
    for file_pair in file_pairs:
        data1 = read_data(os.path.join(folder_path, file_pair[0]))
        data2 = read_data(os.path.join(folder_path, file_pair[1]))
        js_divergence = calculate_js_divergence(data1, data2)  
        output_lines.append(format_output(file_pair, js_divergence))

    output_path = os.path.join(config_output_path, 'Structural-properties-metrics', '3D', 'JSD-of-bond-angles')
    os.makedirs(output_path, exist_ok=True)

    output_file_path = os.path.join(output_path, 'JS-divergence-result.txt')
    with open(output_file_path, 'w') as output_file:
        for line in output_lines:
            output_file.writelines(line)

