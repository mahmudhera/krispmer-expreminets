import pandas as pd
from Bio import SeqIO

summarized_filename = 'mart_export_present_v_absent.txt'
exons_filename = '/home/atif/human_assemblies_kmer_count/exons.fa'

if __name__ == '__main__':
    df_summary = pd.read_csv(summarized_filename, delimiter='\t', header=None)
    print(df_summary)
    transcript_ids = df_summary[3].tolist()

    target_counter = 1

    with open(exons_filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for transcript_id in transcript_ids[:30]:
                if transcript_id in record.id:
                    target_name = f'target_{target_counter}_{transcript_id}.fasta'
                    target_counter += 1
                    f = open(target_name, 'w')
                    f.write('> ' + record.id + '_' + transcript_id + '\n')
                    f.write(str(record.seq) + '\n')
                    f.close()
