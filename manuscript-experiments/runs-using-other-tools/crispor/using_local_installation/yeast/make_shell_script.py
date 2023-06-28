
from os import listdir
from os.path import *

script_file = open('script_for_yeast_crispor.sh', 'w')

targets = [f for f in listdir('./yeast_targets') if isfile(join('./yeast_targets', f))]
print(targets[:10])
for target in targets:
    cmd = f'python /var/www/html/crispor.py sacCer3 yeast_targets/{target} results/{target}.tsv -o results/{target}.tsv'
    script_file.write(cmd + '\n')

script_file.close()
