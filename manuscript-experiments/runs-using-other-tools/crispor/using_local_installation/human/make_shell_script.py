
from os import listdir
from os.path import *

script_file = open('script_chm13_crispor_present_vs_absent.sh', 'w')
target_dir_name = 'chm13_targets_multiple_vs_single'
results_dir_name = 'chm13_results_multiple_vs_single'

targets = [f for f in listdir(target_dir_name) if isfile(join(target_dir_name, f))]
print(targets[:10])
for target in targets:
    cmd = f'python /var/www/html/crispor.py hg38 {target_dir_name}/{target} {results_dir_name}/{target}_guides.tsv'
    script_file.write(cmd + '\n')
script_file.close()
