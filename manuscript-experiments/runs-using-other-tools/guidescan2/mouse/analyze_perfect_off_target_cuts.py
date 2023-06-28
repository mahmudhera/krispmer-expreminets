import dna_jellyfish as jellyfish
import pandas as pd
from os import listdir
from os.path import isfile, join
import subprocess

#jf_file_kr = '23_mer_counts.jf'
#jf_file_gs = '20_mer_counts.jf'

#directory_kr = 'krispmer_targets'
#directory_gs = 'gs_out'
#directory_targets = 'inputs'

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

def generate_gs_out_filename(tgt_name):
    # given a target filename, generate guidescan output filename
    return join(directory_gs, 'gs_out_'+tgt_name)

if __name__ == '__main__':
    mypath = directory_targets
    targets = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]
    print(len(targets))
    print(targets[0])

    ot_count_kr_per_tgt = []
    ot_count_kr_per_grna = []
    ot_count_gs_per_tgt = []
    ot_count_gs_per_grna = []

    qf_genome_gs = jellyfish.QueryMerFile(jf_file_gs)
    qf_genome_kr = jellyfish.QueryMerFile(jf_file_kr)

    for target_file in targets:
        try:
            cmd = 'jellyfish count -m 23 -s 500M -C -o f_23.jf ' + target_file
            args = cmd.split(' ')
            subprocess.call(args)

            only_fname = target_file.split('/')[-1]

            kr_file_with_path = generate_krispmer_out_filename(only_fname)

            df = pd.read_csv(kr_file_with_path)
            df = df[ df['inverse_specificity'] <= cut_off_score ]
            kr_grnas = df['tgt_in_plus'].tolist()

            qf_target = jellyfish.QueryMerFile('f_23.jf')

            total_ot_for_target = 0
            for grna in kr_grnas:
                mer1 = jellyfish.MerDNA(grna)
                mer2 = jellyfish.MerDNA(reverse_complement(grna))
                count_in_target = qf_target[mer1] + qf_target[mer2]
                count_in_genome = qf_genome_kr[mer1] + qf_genome_kr[mer2]
                ot_count = max(0, count_in_genome - count_in_target)

                ot_count_kr_per_grna.append(ot_count)
                total_ot_for_target += ot_count
            ot_count_kr_per_tgt.append(total_ot_for_target)

            print(ot_count_kr_per_tgt)
            print(ot_count_kr_per_grna)

            gs_file_with_path = generate_gs_out_filename(only_fname)

            f = open(gs_file_with_path, 'r')
            all_gs_grnas = f.readlines()
            f.close()

            cmd = 'jellyfish count -m 20 -s 500M -C -o f_20.jf ' + target_file
            args = cmd.split(' ')
            subprocess.call(args)

            qf_target = jellyfish.QueryMerFile('f_20.jf')

            total_ot_for_target = 0
            for grna in all_gs_grnas:
                grna = grna.strip()
                if len(grna) < 23:
                    continue
                if 'CCN' in grna:
                    grna = grna[3:]
                else:
                    grna = grna[:-3]
                mer1 = jellyfish.MerDNA(grna)
                mer2 = jellyfish.MerDNA(reverse_complement(grna))
                count_in_target = qf_target[mer1] + qf_target[mer2]
                count_in_genome = qf_genome_gs[mer1] + qf_genome_gs[mer2]
                ot_count = max(0, count_in_genome - count_in_target)

                ot_count_gs_per_grna.append(ot_count)
                total_ot_for_target += ot_count
            ot_count_gs_per_tgt.append(total_ot_for_target)

            print(ot_count_gs_per_tgt)
            print(ot_count_gs_per_grna)
        except:
            continue

    f = open('perfect_ot_counts_krispmer_per_grna.csv', 'w')
    for x in ot_count_kr_per_grna:
        f.write(str(x))
        f.write('\n')
    f.close()

    f = open('perfect_ot_counts_krispmer_per_target.csv', 'w')
    for x in ot_count_kr_per_tgt:
        f.write(str(x))
        f.write('\n')
    f.close()

    f = open('perfect_ot_counts_guidescan_per_grna.csv', 'w')
    for x in ot_count_gs_per_grna:
        f.write(str(x))
        f.write('\n')
    f.close()

    f = open('perfect_ot_counts_guidecan_per_target.csv', 'w')
    for x in ot_count_gs_per_tgt:
        f.write(str(x))
        f.write('\n')
    f.close()