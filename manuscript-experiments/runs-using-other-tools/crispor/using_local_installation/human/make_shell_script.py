
from os import listdir
from os.path import *

script_file = open('script_chm13_crispor_present_vs_absent.sh', 'w')
target_dir_name = 'chm13_targets_present_vs_absent'
results_dir_name = 'chm13_results_present_vs_absent'

targets = [f for f in listdir('./human_targets') if isfile(join('./human_targets', f))]
print(targets[:10])
for target in targets:
    cmd = f'python /var/www/html/crispor.py hg38 {target_dir_name}/{target} {results_dir_name}/{target}_guides.tsv'
    script_file.write(cmd + '\n')
script_file.close()
