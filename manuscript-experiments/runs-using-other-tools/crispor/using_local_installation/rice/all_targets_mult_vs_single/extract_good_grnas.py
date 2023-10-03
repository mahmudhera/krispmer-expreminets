import pandas as pd
import subprocess

df = pd.read_csv('multiple_v_single.txt', header=None, delimiter='\t')

i = 1
total = len(df[0].tolist())
for kmer in tqdm(df[0].tolist()):
    print(f'Handling {i}/{total}:')
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
