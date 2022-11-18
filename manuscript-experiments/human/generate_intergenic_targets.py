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

    df = pd.read_csv(gbff_filename, comment='#', delimiter='\t', header=None)
    df_intergenes = df[ df[2].str.contains('exon') ]

    bash_file = open('script_for_intergenic_targets.sh', 'w')

    contig_names = df_intergenes[0].tolist()
    start_positions = df_intergenes[3].tolist()
    end_positions = df_intergenes[4].tolist()

    all_genes = list( zip(contig_names, start_positions, end_positions) )
    random.shuffle(all_genes)

    print('indexing genome file...')
    record_dict = SeqIO.index(fasta_filename, "fasta")
    print('indexing complete!')
    all_keys = list(record_dict.keys())

    num_targets_generated = 0
    for (contig_id, start_pos, end_pos) in all_genes:
        if contig_id not in all_keys:
            continue
        start = int(start_pos)
        end = int(end_pos)
        if end < 150 or end-150 <= start:
            continue
        target_start = random.randint(start, end-150)
        seq_record = record_dict[contig_id]
        target_string = str(seq_record.seq[target_start : target_start+150])

        chr_name = contig_id
        fname = 'target_intergene_' + chr_name + "_" + str(target_start) + ':' + str(target_start+150) + '.fasta'
        f = open(fname, 'w')
        f.write('> ' + fname + '\n')
        f.write(target_string)
        f.close()
        bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n -J mer_counts.jf -H k_spectrum_histo\n')

        num_targets_generated += 1
        print('Num of targets generated: ' + str(num_targets_generated))

        if num_targets_generated == num_targets:
            break

    bash_file.close()