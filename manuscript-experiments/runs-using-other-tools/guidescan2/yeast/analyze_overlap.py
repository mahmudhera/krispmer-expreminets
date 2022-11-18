import pandas as pd
from os import listdir
from os.path import isfile, join

targets_dir_name = 'inputs'
gs_out_dir_name = 'gs_out'
kr_out_dir_name = 'krispmer_targets'
cut_off_score = 1.01

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
    print(df.sample(10))
    df = df[ df['inverse_specificity'] <= cut_off_score ]
    print(df.sample(10))
    all_krispmer_targets = df['tgt_in_plus'].tolist() + df['tgt_in_minus'].tolist()


    return None


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
    df = DataFrame(summary, columns=['target_filename', 'only_in_kr', 'in_both', 'only_in_gs'])
    df.to_csv('overlap_summary.csv')
