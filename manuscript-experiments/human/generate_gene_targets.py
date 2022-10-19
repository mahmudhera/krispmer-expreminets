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

    print('reading the gbff file...')
    gbff_filename = 'GCF_000001405.39_GRCh38.p13_genomic.gff'
    parsed_gbff = SeqIO.parse(gbff_filename, 'genbank')
    gbff_records = []
    for rec in parsed_gbff:
        print(rec)
        break
    # read the gbff file
    # read all records into dataframe

    #df.columns = ['acc_id', 'mtype', 'ftype', 'msymbol', 'mname', 'chr', 'start_coord', 'end_coord', 'strand', 'gbuild', 'pcol', 'pdisp']

    #print('Unique chromosomes..')
    #tmp_list = list(set( df['chr'].tolist() ))
    #tmp_list.sort()
    #print( tmp_list )

    #print('Unique feature types...')
    #print( set( df['ftype'].tolist() ) )

    '''
    df_protein_coding_genes = df[ df['ftype'] == 'protein coding gene' ]
    bash_file = open('script_for_protein_coding_targets.sh', 'w')

    print('All seq IDs...')
    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        print(seq_record.id)

    chr_names = df_protein_coding_genes['chr'].tolist()
    start_pos = df_protein_coding_genes['start_coord'].tolist()
    end_pos = df_protein_coding_genes['end_coord'].tolist()

    list_positions = list(zip(chr_names, start_pos, end_pos))
    random.shuffle(list_positions)
    target_id = 1
    count = 0
    for chr_name, start, end in list_positions:
        start = int(start)
        end = int(end)
        if end < 150 or end-150 <= start:
            continue
        target_start = random.randint(start, end-150)
        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == chr_to_contig_id[chr_name]:
                fname = 'target_gene_' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '.fasta'
                f = open(fname, 'w')
                f.write('> ' + 'target_gene_' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '\n')
                f.write(str(seq_record.seq[target_start:target_start+150]))
                f.close()
                target_id = target_id + 1
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n -J mer_counts.jf -H k_spectrum_histo\n')
        count += 1
        if count >= num_targets:
            break

    bash_file.close()
    '''
