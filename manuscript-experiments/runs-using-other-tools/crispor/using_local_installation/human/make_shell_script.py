
from os import listdir
from os.path import *

script_file = open('script_for_human_crispor.sh', 'w')

targets = [f for f in listdir('./human_targets') if isfile(join('./human_targets', f))]
print(targets[:10])
for target in targets:
    cmd = f'python /var/www/html/crispor.py hg38 human_targets/{target} results/{target}.tsv -o results/{target}.tsv'
    script_file.write(cmd + '\n')

script_file.close()
