import argparse
import os
from os import listdir
from os.path import isfile, join
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description="This script will run guidescan on the targets in input directory.\nThe genome needs to be in this directory.\nThe index computed by guidescan neeeds to in this directory\nThe resulting gRNAs will be in results directory.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("guidescan_index_name", type=str, help="The name of the guidescan index, which needs to be in this directory")
    parser.add_argument("input_dir", type=str, help="Directory where all targets are")
    parser.add_argument("output_dir", type=str, help="Directory where all results will be. This must already exist before running this script.")
    args = parser.parse_args()
    return args.guidescan_index_name, args.input_dir, args.output_dir

if __name__ == "__main__":
    index_name, input_dir, output_dir = parse_args()

    # get all targets as list
    mypath = input_dir
    targets = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]

    for target_file in targets:
        # first creake kmers file
        f = open('tmp_kmer_file', 'w')
        cmd = 'python generate_kmers.py --kmer-length 20 ' + target_file + ' --min-chr-length 1'
        args = cmd.split(' ')
        subprocess.call(args, stdout=f)
        f.close()

        # then call guidescan enumerate
        out_filename = join(output_dir, 'gs_out_' + target_file.split('/')[-1].split('.fasta')[0])
        print(out_filename)
        cmd = 'guidescan enumerate -m 3 -f tmp_kmer_file ' + index_name + ' --output ' + out_filename
        args = cmd.split(' ')
        subprocess.call(args)
