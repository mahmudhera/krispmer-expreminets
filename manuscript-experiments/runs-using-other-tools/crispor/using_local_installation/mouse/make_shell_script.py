
from os import listdir
from os.path import *

script_file = open('script_for_mouse_crispor.sh', 'w')

targets = [f for f in listdir('./mouse_targets') if isfile(join('./mouse_targets', f))]
print(targets[:10])
for target in targets:
    cmd = f'python /var/www/html/crispor.py mm9 mouse_targets/{target} results/{target}.tsv -o results/{target}.tsv'
    script_file.write(cmd + '\n')

script_file.close()
