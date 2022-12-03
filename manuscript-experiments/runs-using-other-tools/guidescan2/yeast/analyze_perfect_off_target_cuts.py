import jellyfish
import pandas as pd
from os import listdir
from os.path import isfile, join
import subprocess

jf_file_kr = '23_mer_counts.jf'
jf_file_gs = '20_mer_counts.jf'

directory_kr = 'krispmer_targets'
directory_gs = 'gs_out'
directory_targets = 'inputs'

cut_off_score = 1.2

def complement(seq):
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)

def reverse_complement(s):
    return complement(s[::-1])

def generate_krispmer_out_filename(tgt_name):
    # given a target filename, generate krispmer output filename
    return join(directory_kr, 'scores_'+tgt_name)

if __name__ == '__main__':
    mypath = directory_targets
    targets = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]
    print(len(targets))
    print(targets[0])

    ot_count_kr_per_tgt = []
    ot_count_kr_per_grna = []
    ot_count_gs_per_tgt = []
    ot_count_gs_per_grna = []

    for target_file in targets:
        cmd = 'jellyfish count -m 23 -s 500M -C -o f_23.jf ' + target_file
        args = cmd.split(' ')
        subprocess.call(args)

        only_fname = target_file.split('/')[-1]

        kr_file_with_path = generate_krispmer_out_filename(only_fname)

        df = pd.read_csv(kr_file_with_path)
        df = df[ df['inverse_specificity'] <= cut_off_score ]
        kr_grnas = df['tgt_in_plus'].tolist()

        qf_target = jellyfish.QueryMerFile('f_23.jf')
        qf_genome = qf = jellyfish.QueryMerFile(jf_file_kr)

        total_ot_for_target = 0
        for grna in kr_grnas:
            mer1 = jellyfish.MerDNA(grna)
            mer2 = jellyfish.MerDNA(reverse_complement(grna))
            count_in_target = qf_target[mer1] + qf_target[mer2]
            count_in_genome = qf_genome[mer1] + qf_genome[mer2]
            ot_count = max(0, count_in_genome - count_in_target)

            ot_count_kr_per_grna.append(ot_count)
            total_ot_for_target += ot_count
        ot_count_kr_per_tgt.append(total_ot_for_target)

        print(ot_count_kr_per_tgt)
        print(ot_count_kr_per_grna)

        exit(-1)

    # for each input file e:
        # generate 23 count file -> f_23
        # locate the kr file -> kr
        # for each in kr file:
            # filter using threshold
            # find count in genome (23)
            # find count in f_23
            # find ot_count
            # record count accordingly
        # generate 20 count file -> f_20
        # locate the gs file -> kr
        # for each in gs file:
            # find count in genome (20)
            # find count in f_23
            # find ot_count
            # record count accordingly
    # report all counts
