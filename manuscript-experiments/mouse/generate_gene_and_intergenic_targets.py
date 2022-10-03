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


chr_to_contig_id = {
    '1' : 'NC_000067.7',
    '2' : 'NC_000068.8',
    '3' : 'NC_000069.7',
    '4' : 'NC_000070.7',
    '5' : 'NC_000071.7',
    '6' : 'NC_000072.7',
    '7' : 'NC_000073.7',
    '8' : 'NC_000074.7',
    '9' : 'NC_000075.7',
    '10' : 'NC_000076.7',
    '11' : 'NC_000077.7',
    '12' : 'NC_000078.7',
    '13' : 'NC_000079.7',
    '14' : 'NC_000080.7',
    '15' : 'NC_000081.7',
    '16' : 'NC_000082.7',
    '17' : 'NC_000083.7',
    '18' : 'NC_000084.7',
    '19' : 'NC_000085.7',
    'X' : 'NC_000086.8',
    'Y' : 'NC_000087.8'
}


if __name__ == "__main__":
    random.seed(random_seed)
    df = pd.read_csv('MGI_MRK_Coord.tsv', delimiter='\t', index_col=False, dtype=str)

    # the following are the headers
    # ['1. MGI Marker Accession ID', '2. Marker Type', '3. Feature Type','4. Marker Symbol', '5. Marker Name', '6. Chromosome','7. Start Coordinate', '8. End Coordinate', '9. Strand','10. Genome Build', '11. Provider Collection', '12. Provider Display']
    # lets make them easier strings
    df.columns = ['acc_id', 'mtype', 'ftype', 'msymbol', 'mname', 'chr', 'start_coord', 'end_coord', 'strand', 'gbuild', 'pcol', 'pdisp']

    print('Unique chromosomes..')
    tmp_list = list(set( df['chr'].tolist() ))
    tmp_list.sort()
    print( tmp_list )

    print('Unique feature types...')
    print( set( df['ftype'].tolist() ) )

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
        if end < 150 or end-150 <= start:
            continue
        target_start = random.randint(start, end-150)
        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == dic[chr_name]:
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

    # intergenic regions
    df_intergenic = df[(df['chr']!='protein coding gene')]
    bash_file = open('script_for_intergenic_targets.sh', 'w')

    chr_names = df_intergenic['chr'].tolist()
    start_pos = df_intergenic['start_coods'].tolist()
    end_pos = df_intergenic['end_coord'].tolist()

    list_positions = list(zip(chr_names, start_pos, end_pos))
    random.shuffle(list_positions)
    count = 0
    target_id = 0
    for chr_name, start, end in list_positions:
        if end < 150 or end-150 <= start:
            continue
        target_start = random.randint(start, end-150)
        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == dic[chr_name]:
                fname = 'target_intergenic_' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '.fasta'
                f = open(fname, 'w')
                f.write('> ' + 'target_intergenic_' + str(target_id) + '_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '\n')
                f.write(str(seq_record.seq[target_start:target_start+150]))
                f.close()
                target_id = target_id + 1
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n -J mer_counts.jf -H k_spectrum_histo\n')
        count += 1
        if count >= num_targets:
            break

    bash_file.close()
