import pandas as pd
from Bio import SeqIO
import random

random.seed(0)

candidate_genes_list_filename = '/home/atif/oryza-sativa-genomes-and-counts/kmer_counts/genes_multiple_v_single_analysis_all.txt'
sequences_filename = '/home/atif/oryza-sativa-genomes-and-counts/kmer_counts/cds_primary.fa'

if __name__ == '__main__':
    # read all candidate genes list
    f = open(candidate_genes_list_filename, 'r')
    lines = f.readlines()
    f.close()
    gene_names_list = [ line.strip() for line in lines ]
    random.shuffle(gene_names_list)
    gene_names_list = gene_names_list[:100]

    target_counter = 1

    with open(sequences_filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for gene_name in gene_names_list:
                if gene_name in record.id:
                    target_name = f'target_{target_counter}_{gene_name}.fasta'
                    target_counter += 1
                    f = open(target_name, 'w')
                    f.write('>' + target_name + '\n')
                    f.write(str(record.seq) + '\n')
                    f.close()
