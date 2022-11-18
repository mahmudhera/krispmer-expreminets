import pandas as pd
from os import listdir
from os.path import isfile, join

targets_dir_name = 'inputs'
gs_out_dir_name = 'gs_out'
kr_out_dir_name = 'krispmer_targets'
cut_off_score = 1.2

def generate_gs_out_filename(tgt_name):
    # given a target filename, generate guidescan output filename
    return join(gs_out_dir_name, 'gs_out_'+tgt_name)


def generate_krispmer_out_filename(tgt_name):
    # given a target filename, generate krispmer output filename
    return join(kr_out_dir_name, 'scores_'+tgt_name)


def find_overlap_given_target_filename(target_filename):
    krispmer_filename = generate_krispmer_out_filename(target_filename)
    guidescan_filename = generate_gs_out_filename(target_filename)

    df = pd.read_csv(krispmer_filename)
    df = df[ df['inverse_specificity'] <= cut_off_score ]
    all_krispmer_grnas = df['tgt_in_plus'].tolist() + df['tgt_in_minus'].tolist()

    f = open(guidescan_filename, 'r')
    all_gs_grnas = f.readlines()
    f.close()

    only_in_kr, in_both, only_in_gs = 0,0,0
    for potential_grna in all_gs_grnas:
        if len(potential_grna) < 20:
            continue

        only_in_gs += 1
        gs_grna = potential_grna.upper()
        if 'CCN' in gs_grna:
            gs_grna = gs_grna[3:]
        else:
            gs_grna = gs_grna[:-3]

        found = False
        for kr_grna in all_krispmer_grnas:
            if gs_grna in kr_grna:
                found = True
                break

        if found:
            in_both += 1

    only_in_gs -= in_both
    only_in_kr = len(all_krispmer_grnas)/2 - in_both

    return only_in_kr, in_both, only_in_gs


def find_all_targets():
    mypath = targets_dir_name
    targets = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]
    return targets


if __name__ == '__main__':
    target_filenames_list = find_all_targets()
    summary = []
    for target_filename in target_filenames_list:
        target_name = str(target_filename.split('/')[-1])
        print(target_name)
        print( generate_gs_out_filename(target_name) )
        print( generate_krispmer_out_filename(target_name) )
        only_krispmer, common, only_guidescan = find_overlap_given_target_filename(target_name)
        summary.append( (target_filename, only_krispmer, common, only_guidescan) )
    df = pd.DataFrame(summary, columns=['target_filename', 'only_in_kr', 'in_both', 'only_in_gs'])
    df.to_csv('overlap_summary.csv')
