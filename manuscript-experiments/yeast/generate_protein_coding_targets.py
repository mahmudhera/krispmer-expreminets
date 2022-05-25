import numpy as np
import pandas as pd
import random
from Bio import SeqIO

num_targets = 50
target_length = 150
random_seed = 1
fasta_filename = 'genome.fasta'
hd = 3

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
    for chr_name, start, end in list_positions[:num_targets]:
        target_start = random.randint(start, end-150)
        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == dic[chr_name]:
                print(chr_name, seq_record.id, target_start, target_start+150)



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
    for chr_name, start, end in list_positions[:num_targets]:
        break
        if end < 150 or end-150 <= start:
            continue
        target_start = random.randint(start, end-150)
        print(chr_name, target_start, target_start+150)
