import pandas as pd
from os import listdir

targets_dir_name = 'inputs'
gs_out_dir_name = 'gs_out'
kr_out_dir_name = 'krispmer_targets'

def generate_gs_out_filename():
    # given a target filename, generate guidescan output filename
    return None


def generate_krispmer_out_filename():
    # given a target filename, generate krispmer output filename
    return None


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
    print(target_filenames_list[:5])
    print(target_filenames_list[-5:])
    exit(-1)
    summary = []
    for target_filename in target_filenames_list:
        only_krispmer, common, only_guidescan = find_overlap_given_target_filename(target_filename)
        summary.append( (target_filename, only_krispmer, common, only_guidescan) )
    df = DataFrame(summary, columns=['target_filename', 'only_in_kr', 'in_both', 'only_in_gs'])
    df.to_csv('overlap_summary.csv')
