import pandas as pd
from Bio import SeqIO
import random

random.seed(0)

exons_filename = '/home/atif/human_assemblies_kmer_count/exons.fa'
handpicked_transcript_ids = ['ENST00000375050.6', 'ENST00000288050.9', 'ENST00000261405.10', 'ENST00000351205.8', 'ENST00000380742.8', 'ENST00000374312.5']

if __name__ == '__main__':
    target_counter = 1

    with open(exons_filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for transcript_id in handpicked_transcript_ids:
                if transcript_id in record.id:
                    target_name = f'target_{target_counter}_{transcript_id}.fasta'
                    target_counter += 1
                    f = open(target_name, 'w')
                    f.write('>' + target_name + '\n')
                    f.write(str(record.seq) + '\n')
                    f.close()

    with open(exons_filename) as handle:
        all_record_ids = []
        for record in SeqIO.parse(handle, "fasta"):
            all_record_ids.append(str(record.id))
        random.shuffle(all_record_ids)
        selected_record_ids = all_record_ids[:100]

        print(selected_record_ids)
    with open(exons_filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if str(record.id) in selected_record_ids:
                target_name = f'target_{target_counter}_{record.id}.fasta'
                target_counter += 1
                f = open(target_name, 'w')
                f.write('>' + target_name + '\n')
                f.write(str(record.seq) + '\n')
                f.close()
