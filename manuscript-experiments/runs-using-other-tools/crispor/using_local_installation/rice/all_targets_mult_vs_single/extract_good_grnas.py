import pandas as pd
import subprocess
import os
import random

df = pd.read_csv('multiple_v_single.txt', header=None, delimiter='\t')

f_out = open('grnas_and_spec_scores', 'w')
f_out.write('grna \t spec_score\n')

i = 1
total = len(df[0].tolist())
for kmer in random.shuffle(df[0].tolist()):
    print('Handling '+str(i)+'/'+str(total)+':')
    i += 1
    f = open('tmp.fasta', 'w')
    f.write('>target\n')
    f.write(kmer+'\n')
    f.close()
    cmd = 'python /var/www/html/crispor.py pz9Osativa tmp.fasta tmp.tsv'
    subprocess.call(cmd.split(' '), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    f = open('tmp.tsv', 'r')
    line = f.readlines()[1]
    f.close()
    # extract the specificity score
    spec_score = float( line.split('\t')[3] )
    print(kmer, spec_score)
    f_out.write(kmer + '\t' + str(spec_score) + '\n')
    f_out.flush()

f_out.close()
