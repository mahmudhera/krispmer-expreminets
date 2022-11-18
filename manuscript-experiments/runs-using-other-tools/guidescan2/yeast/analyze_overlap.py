import pandas as pd
from os import listdir
from os.path import isfile, join

targets_dir_name = 'inputs'
gs_out_dir_name = 'gs_out'
kr_out_dir_name = 'krispmer_targets'

def generate_gs_out_filename(tgt_name):
    # given a target filename, generate guidescan output filename
    return join(gs_out_dir_name, 'gs_out_'+tgt_name)


def generate_krispmer_out_filename(tgt_name):
    # given a target filename, generate krispmer output filename
    return join(kr_out_dir_name, 'scores_'+tgt_name)


def find_overlap_given_target_filename(target_filename):
    # given a target file, locate the krispmer and the guidescan output file
    # then, return overlap
    return None


def find_all_targets():
    mypath = targets_dir_name
    targets = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]
    return targets


if __name__ == '__main__':
    target_filenames_list = find_all_targets()
    summary = []
    for target_filename in target_filenames_list:
        target_name = target_filename.split('/')[-1]
        print( generate_gs_out_filename(target_name) )
        print( generate_krispmer_out_filename(target_name) )
        only_krispmer, common, only_guidescan = find_overlap_given_target_filename(target_filename)
        summary.append( (target_filename, only_krispmer, common, only_guidescan) )
    df = DataFrame(summary, columns=['target_filename', 'only_in_kr', 'in_both', 'only_in_gs'])
    df.to_csv('overlap_summary.csv')
