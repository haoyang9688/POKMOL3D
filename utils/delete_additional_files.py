import os
import re

def delete_additional_files(config):
    directory = config['Code_path']
    
    pattern_starts_with_server = re.compile(r'.*(localhost-|server-).*')
    pattern_ends_with_log = re.compile(r'\d+\.log$')
    pattern_ends_with_others = re.compile(r'.*\.(log|sdf\.log|maegz|sdf|txt)$')
    
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
      
        if os.path.isfile(file_path):
            if (pattern_starts_with_server.match(file_name) or 
                pattern_ends_with_log.match(file_name) or 
                pattern_ends_with_others.match(file_name)):
                #print(f"Deleting file: {file_path}")
                os.remove(file_path)




