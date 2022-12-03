import jellyfish
import pandas as pd
from os import listdir
from os.path import isfile, join

jf_file_kr = '23_mer_counts.jf'
jf_file_gs = '20_mer_counts.jf'

directory_kr = 'krispmer_targets'
directory_gs = 'gs_out'
directory_targets = 'inputs'

if __name__ == '__main__':
    mypath = directory_targets
    targets = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]
    print(len(targets))
    print(targets[0])
    # ot_count_kr_per_tgt = []
    # ot_count_kr_per_grna = []
    # ot_count_gs_per_tgt = []
    # ot_count_gs_per_grna = []
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
