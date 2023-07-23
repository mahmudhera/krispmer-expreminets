import pandas as pd

summarized_filename = '/home/atif/human_assemblies_kmer_count/multiple_v_single_candidates.txt'
exons_filename = '/home/atif/human_assemblies_kmer_count/exons.fa'

if __name__ == '__main__':
    df_summary = pd.read_csv(summarized_filename, delimiter='\t', header=None)
    print(df_summary)
