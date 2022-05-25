import numpy as np
import pandas as pd
import random
from Bio import SeqIO

num_targets = 50
target_length = 150
random_seed = 1
fasta_filename = 'genome.fasta'
hd = 3
reads_filename = 'reads.fastq'

dic = {}
dic['chrI'] = 'NC_001133.9'
dic['chrII'] = 'NC_001134.8'
dic['chrIII'] = 'NC_001135.5'
dic['chrIV'] = 'NC_001136.10'
dic['chrV'] = 'NC_001137.3'
dic['chrVI'] = 'NC_001138.5'
dic['chrVII'] = 'NC_001139.9'
dic['chrVIII'] = 'NC_001140.6'
dic['chrIX'] = 'NC_001141.2'
dic['chrX'] = 'NC_001142.9'
dic['chrXI'] = 'NC_001143.9'
dic['chrXII'] = 'NC_001144.5'
dic['chrXIII'] = 'NC_001145.3'
dic['chrXIV'] = 'NC_001146.8'
dic['chrXV'] = 'NC_001147.6'
dic['chrXVI'] = 'NC_001148.4'
dic['chrmt'] = 'NC_001224.1'

if __name__ == "__main__":
    random.seed(random_seed)

    bash_file = open('script.sh', 'w')

    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        print(seq_record.id)

    # non protein coding genes
    df = pd.read_csv('AGR_Yeast_Genes.tsv', delimiter='\t', header=None)
    df.columns = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
    df = df[df['j']==1]
    df = df[df['l']=='ORF']

    chr_names = df['g'].tolist()
    start_pos = df['h'].tolist()
    end_pos = df['i'].tolist()

    list_positions = list(zip(chr_names, start_pos, end_pos))
    random.shuffle(list_positions)
    target_id = 1
    count = 0
    for chr_name, start, end in list_positions:
        if end < 150 or end-150 <= start:
            continue
        target_start = random.randint(start, end-150)
        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == dic[chr_name]:
                fname = 'target' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '.fasta'
                f = open(fname, 'w')
                f.write('> ' + 'target' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '\n')
                f.write(str(seq_record.seq[target_start:target_start+150]))
                f.close()
                target_id = target_id + 1
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n\n')
        count += 1
        if count == num_targets:
            break


    # protein coding genes
    df = pd.read_csv('AGR_Yeast_Genes.tsv', delimiter='\t', header=None)
    df.columns = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
    df = df[df['j']==1]
    df = df[(df['l']=='transposable element gene') | (df['l']=='snoRNA gene') | (df['l']=='snRNA gene') | (df['l']=='rRNA gene') | (df['l']=='ncRNA gene') | (df['l']=='tRNA gene') | (df['l']=='telomerase RNA gene')]

    chr_names = df['g'].tolist()
    start_pos = df['h'].tolist()
    end_pos = df['i'].tolist()

    list_positions = list(zip(chr_names, start_pos, end_pos))
    random.shuffle(list_positions)
    count = 0
    for chr_name, start, end in list_positions:
        count += 1
        if end < 150 or end-150 <= start:
            continue
        target_start = random.randint(start, end-150)
        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == dic[chr_name]:
                fname = 'target' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '.fasta'
                f = open(fname, 'w')
                f.write('> ' + 'target' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '\n')
                f.write(str(seq_record.seq[target_start:target_start+150]))
                f.close()
                target_id = target_id + 1
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n\n')
        count += 1
        if count == num_targets:
            break

    bash_file.close()
