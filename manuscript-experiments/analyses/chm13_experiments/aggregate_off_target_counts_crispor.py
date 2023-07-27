import argparse
import os
from itertools import chain, combinations, product
import pandas as pd

specificity_cutoff = 60

def complement(seq):
    """
    generates the complement sequence, e.g.: ACGT-->TCGA
    :param seq: dna sequence
    :return: complement
    """
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(str(seq))
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)

def reverse_complement(s):
    return complement(s[::-1])

def get_all_files(directory):
    files = os.listdir(directory)
    for file in files:
        if os.path.isfile(os.path.join(directory, file)):
            yield os.path.join(directory, file)

def generate_scores_filename(target_filename_with_path):
    splitted = target_filename_with_path.split('/')
    return '/'.join(splitted[:-2] + ['crispor_guides', splitted[-1]+'_guides.tsv'])

def generate_off_target_filename(target_filename_with_path):
    splitted = target_filename_with_path.split('/')
    return '/'.join(splitted[:-2] + ['crispor_grnas', splitted[-1]+'_off_targets.tsv'])

def main():
    # take command line arguments
    parser = argparse.ArgumentParser(description='Process guidescan2 results.')
    parser.add_argument('target_dir_name', help='Full path to the directory containing all target sequences')
    parser.add_argument('grna_dir_name', help='Full path to the directory containing all grnas')

    args = parser.parse_args()

    target_dir_name = args.target_dir_name
    grna_dir_name = args.grna_dir_name

    target_files_and_grna_files = []
    for target_file in get_all_files(target_dir_name):
        target_files_and_grna_files.append((target_file, generate_scores_filename(target_file), generate_off_target_filename(target_file)))
    # list all target and grna files

    for target_file, grna_file, off_target_file in target_files_and_grna_files:
        df_grnas = pd.read_csv(grna_file, delimiter='\t')
        df_off_targets = pd.read_csv(off_target_file, delimiter='\t')

        grnas = df_grnas['guideId'].tolist()
        off_target_counts = df_grnas['mitSpecScore'].tolist()
        specificity_scores = df_grnas['targetSeq'].tolist()

        for grna, off_target_count, specificity_score in list(zip(grnas, off_target_counts, specificity_scores)):
            if str(specificity_score) == 'None':
                continue
            if float(specificity_score) < specificity_cutoff:
                continue
            df_this_grna = df_off_targets[ df_off_targets['guideSeq']==grna ]
            num_mismatches_list = df_this_grna['mismatchCount'].tolist()
            a = num_mismatches_list.count(0)
            b = num_mismatches_list.count(1)
            c = int(off_target_count) - (a+b)

            ot_count_0_mismatch, ot_count_1_mismatch, ot_count_2_mismatch = a,b,c
            print(str(target_file).split('/')[-1], grna, ot_count_0_mismatch, ot_count_1_mismatch, ot_count_2_mismatch)

if __name__ == '__main__':
    main()
