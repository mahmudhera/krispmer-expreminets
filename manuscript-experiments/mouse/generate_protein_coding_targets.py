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

if __name__ == "__main__":
    random.seed(random_seed)
    df = pd.read_csv('MGI_MRK_Coord.rpt', delimiter='\t')
    print(df)
    '''
    bash_file = open('script_protein_coding.sh', 'w')

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
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n -J mer_counts.jf -H k_spectrum_histo\n')
        count += 1
        if count >= num_targets:
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
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n -J mer_counts.jf -H k_spectrum_histo\n')
        count += 1
        if count >= num_targets:
            break

    bash_file.close()
    '''
