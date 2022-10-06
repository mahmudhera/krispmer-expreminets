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
    df = pd.read_csv('mouse_gRNAs_multiple.sam', delimiter='\t', header=None)
    random.seed(random_seed)

    contig_names = df[2].tolist()
    start_pos = df[3].tolist()
    list_positions = list( zip(contig_names, start_pos) )

    # test start
    contig_name_test = contig_names[0]
    start_pos_test = start_pos[0]
    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        if seq_record.id == contig_name_test:
            print( str( seq_record.seq[start_pos_test-10:start_pos_test+30] ) )
    # test end

    random.shuffle(list_positions)
    num_targets_generated = 0

    bash_file = open('script_duplicated_regions.sh', 'w')

    target_count_id = 1
    for position in list_positions:
        contig_name = position[0]
        start_pos = position[1]
        if start_pos < 80:
            continue

        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == contig_name:
                fname = 'target_duplicated_' + str(target_count_id) + '_' + contig_name + "_" + str(start_pos-75) + ':' + str(start_pos+75) + '.fasta'
                f = open(fname, 'w')
                f.write('> ' + 'target_duplicated_' + str(target_count_id) + '_' + contig_name + "_" + str(start_pos-75) + ':' + str(start_pos+75) + '\n')
                target_count_id += 1
                f.write(str(seq_record.seq[start_pos-75:start_pos+75]).upper())
                f.close()
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n -J mer_counts.jf -H k_spectrum_histo\n')

        num_targets_generated+=1
        if num_targets_generated == num_targets:
            break

    bash_file.close()
