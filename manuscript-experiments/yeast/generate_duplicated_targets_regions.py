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
    df = pd.read_csv('yeast_gRNAs_multiple.sam', delimiter='\t', header=None)
    random.seed(random_seed)

    contig_names = df[2].tolist()
    start_pos = df[3].tolist()
    list_positions = list( zip(contig_names, start_pos) )

    random.shuffle(list_positions)
    num_targets_generated = 0

    bash_file = open('script_duplicated.sh', 'w')

    for position in list_positions:
        contig_name = position[0]
        start_pos = position[1]
        if start_pos < 80:
            continue

        for seq_record in SeqIO.parse(fasta_filename, "fasta"):
            if seq_record.id == contig_name:
                fname = 'target_' + chr_name + "_" + str(target_start-75) + ':' + str(target_start+75) + '.fasta'
                f = open(fname, 'w')
                f.write('> ' + 'target_' + chr_name + "_" + str(target_start-75) + ':' + str(target_start+75) + '\n')
                f.write(str(seq_record.seq[target_start-75:target_start+75]))
                f.close()
                bash_file.write('krispmer ' + reads_filename + ' ' + fname + ' ' + 'scores_' + fname + ' ' + str(hd) + ' -n -J mer_counts.jf -H k_spectrum_histo -n\n')

    bash_file.close()
