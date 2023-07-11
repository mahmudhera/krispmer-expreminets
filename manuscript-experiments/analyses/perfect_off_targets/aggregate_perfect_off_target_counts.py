import argparse
import os
import dna_jellyfish as jellyfish

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
        target_file = open(target_file, 'r')
        lines = target_file.readlines()
        target_file.close()
        target_sequence = ''.join( [line.strip().upper() for line in lines[1:]] )
        print(target_sequence)

        for tgt_in_plus, tgt_in_minus in grna_list:
            num_occurrences_in_target = target_sequence.count(tgt_in_plus)+target_sequence.count(tgt_in_minus)
            num_occurrences_in_genome = qf_genome[jellyfish.MerDNA(tgt_in_plus)] + qf_genome[jellyfish.MerDNA(tgt_in_minus)]
            ot_count = max(0, num_occurrences_in_genome - num_occurrences_in_target)
            print(tgt_in_plus, ot_count)
            # <target_filename, grna, ot_count, type> add this to the summary file

if __name__ == '__main__':
    main()
