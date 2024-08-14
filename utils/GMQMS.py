import os

def extract_last_value(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            if lines:
                last_line = lines[-1].strip()
                value = last_line.split(':')[-1].strip()
                return value
            else:
                return 'NA'
    except FileNotFoundError:
        return 'NA'

def process_quality_metrics(base_path, output_file):
    metrics = [
        ('Validity', 'Validity'),
        ('Usability', 'Usability'),
        ('Uniqueness', 'Uniqueness'),
        ('QED', 'QED-mean.txt'),
        ('SA', 'SA-mean.txt'),
        ('Morgan-Interal-diversity', 'morgan-tdiv-mean.txt'),
        ('Scaffold-Interal-diversity', 'scaffold-tdiv-mean.txt')
    ]

    results = ['General molecular quality metrics\tValues\n']

    for metric, file_name in metrics:
        metric_path = os.path.join(base_path, metric)
        if metric in ['QED', 'SA']:
            file_path = os.path.join(metric_path, file_name)
        elif metric in ['Morgan-Interal-diversity', 'Scaffold-Interal-diversity']:
            file_path = os.path.join(metric_path, file_name)
        else:
            # For Validity, Usability, and Uniqueness, find .txt files in the directory
            txt_files = [f for f in os.listdir(metric_path) if f.endswith('.txt')]
            if txt_files:
                file_path = os.path.join(metric_path, txt_files[0])
            else:
                file_path = None

        if file_path:
            value = extract_last_value(file_path)
        else:
            value = 'NA'

        results.append(f"{metric}\t{value}\n")

    with open(output_file, 'w') as out_file:
        out_file.writelines(results)

def calculate_Summary(config):
    config_output_path = config['output']['output_path']
    base_path = os.path.join(config_output_path, 'General-molecular-quality-metrics')
    
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    
    output_file = os.path.join(base_path, 'General-molecular-quality-metrics.txt')
    process_quality_metrics(base_path, output_file)
