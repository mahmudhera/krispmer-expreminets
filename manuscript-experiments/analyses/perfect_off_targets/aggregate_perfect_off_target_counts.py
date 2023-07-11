import argparse
import os
import dna_jellyfish as jellyfish
import trie
from itertools import chain, combinations, product

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

def hamming_circle(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.
    >>> sorted(hamming_circle('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_circle('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_circle('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']
    """
    return_list = []
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            return_list.append(''.join(cousin))
    return return_list

def generate_adjacent_mers(sequence, hamming_distance):
	alphabet = 'AGCT'
	return hamming_circle(sequence, hamming_distance, alphabet)

def get_all_files(directory):
    files = os.listdir(directory)
    for file in files:
        if os.path.isfile(os.path.join(directory, file)):
            yield os.path.join(directory, file)

def generate_scores_filename(target_filename_with_path):
    splitted = target_filename_with_path.split('/')
    return '/'.join(splitted[:-2] + ['krispmer_guides','scores_'+splitted[-1]])

def main():
    # take command line arguments
    parser = argparse.ArgumentParser(description='Process command line arguments.')
    parser.add_argument('target_dir_name', help='Full path to the directory containing all target sequences')
    parser.add_argument('grna_dir_name', help='Full path to the directory containing all grnas')
    parser.add_argument('genome_jf_filename', help='Full path to the jellyfish counted file for the genome')

    args = parser.parse_args()

    target_dir_name = args.target_dir_name
    grna_dir_name = args.grna_dir_name
    genome_jf_filename = args.genome_jf_filename

    # jellyfish query file
    qf_genome = jellyfish.QueryMerFile(genome_jf_filename)

    target_files_and_grna_files = []
    for target_file in get_all_files(target_dir_name):
        target_files_and_grna_files.append((target_file, generate_scores_filename(target_file)))
    # list all target and grna files

    for target_file, grna_file in target_files_and_grna_files:
        # list all grnas
        grna_file = open(grna_file, 'r')
        lines = grna_file.readlines()
        grna_file.close()
        grna_list = []
        for line in lines[1:]:
            if float(line.split(',')[2]) > 1.5:
                continue
            grna_list.append( (line.split(',')[0], line.split(',')[1]) )
        print(grna_list)

        # get the target sequence
        tgt_f = open(target_file, 'r')
        lines = tgt_f.readlines()
        tgt_f.close()
        target_sequence = ''.join( [line.strip().upper() for line in lines[1:]] )
        print(target_sequence)

        for tgt_in_plus, tgt_in_minus in grna_list:
            num_occurrences_in_target = target_sequence.count(tgt_in_plus)+target_sequence.count(tgt_in_minus)
            num_occurrences_in_genome = qf_genome[jellyfish.MerDNA(tgt_in_plus)] + qf_genome[jellyfish.MerDNA(tgt_in_minus)]
            ot_count_0_mismatch = max(0, num_occurrences_in_genome - num_occurrences_in_target)

            sequences_with_one_distance = set(generate_adjacent_mers(tgt_in_plus, 1))
            num_occurrences_in_target = 0
            num_occurrences_in_genome = 0
            for potential_off_target in sequences_with_one_distance:
                num_occurrences_in_target += target_sequence.count(potential_off_target) + target_sequence.count(reverse_complement(potential_off_target))
                num_occurrences_in_genome += qf_genome[jellyfish.MerDNA(potential_off_target)] + qf_genome[jellyfish.MerDNA(reverse_complement(potential_off_target))]
            ot_count_1_mismatch = max(0, num_occurrences_in_genome - num_occurrences_in_target)

            sequences_with_two_distance = set(generate_adjacent_mers(tgt_in_plus, 2))
            num_occurrences_in_target = 0
            num_occurrences_in_genome = 0
            for potential_off_target in sequences_with_one_distance:
                num_occurrences_in_target += target_sequence.count(potential_off_target) + target_sequence.count(reverse_complement(potential_off_target))
                num_occurrences_in_genome += qf_genome[jellyfish.MerDNA(potential_off_target)] + qf_genome[jellyfish.MerDNA(reverse_complement(potential_off_target))]
            ot_count_2_mismatch = max(0, num_occurrences_in_genome - num_occurrences_in_target)

            print(str(target_file).split('/')[-1], tgt_in_plus, ot_count_0_mismatch, ot_count_1_mismatch, ot_count_2_mismatch)
            # <target_filename, grna, ot_count, type> add this to the summary file

if __name__ == '__main__':
    main()
