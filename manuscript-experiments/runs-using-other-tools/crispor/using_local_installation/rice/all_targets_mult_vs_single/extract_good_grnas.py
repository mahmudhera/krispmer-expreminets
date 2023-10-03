import pandas as pd
from tqdm import tqdm
import subprocess

df = pd.read_csv('multiple_v_single.txt', header=None, delimiter='\t')

for kmer in tqdm(df[0].tolist()):
    f = open('tmp.fasta', 'w')
    f.write('>target\n')
    f.write(kmer+'\n')
    f.close()
    cmd = 'python /var/www/html/crispor.py pz9Osativa tmp.fasta tmp.tsv'
    subprocess.run(cmd.split(' '), stdout=open(os.devnull, 'wb'))
    f = open('tmp.tsv', 'r')
    line = f.readlines()[1]
    f.close()
    # extract the specificity score
    spec_score = float( line.split('\t')[3] )
    print(kmer, spec_score)
