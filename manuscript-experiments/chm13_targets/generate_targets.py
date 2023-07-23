import pandas as pd
from Bio import SeqIO

summarized_filename = '/home/atif/human_assemblies_kmer_count/multiple_v_single_candidates.txt'
exons_filename = '/home/atif/human_assemblies_kmer_count/exons.fa'

if __name__ == '__main__':
    df_summary = pd.read_csv(summarized_filename, delimiter='\t', header=None)
    transcript_ids = df_summary[3].tolist()

    target_counter = 1

    with open(exons_filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for transcript_id in transcript_ids:
                if transcript_id in record.id:
                    target_name = f'target_{target_counter}.fasta'
                    target_counter += 1
                    f = open(target_name, 'w')
                    f.write('> ' + record.id + '\n')
                    f.write(record.seq)
                    f.close()
