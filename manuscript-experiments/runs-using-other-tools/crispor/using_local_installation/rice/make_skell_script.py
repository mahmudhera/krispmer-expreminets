from os import listdir
from os.path import *

script_file = open('script_rice.sh', 'w')
target_dir_name = 'multiple_vs_single_targets'
results_dir_name = 'multiple_vs_single_results'

targets = [f for f in listdir(target_dir_name) if isfile(join(target_dir_name, f))]
print(targets[:10])
for target in targets:
    cmd = f'python /var/www/html/crispor.py pz9Osativa {target_dir_name}/{target} {results_dir_name}/{target}_guides.tsv'
    script_file.write(cmd + '\n')
script_file.close()
