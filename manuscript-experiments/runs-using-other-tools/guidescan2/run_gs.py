import argparse
import os
from os import listdir
from os.path import isfile, join

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
    print(targets)

    # for everything in that list
        # first creake kmers file
        # then call guidescan enumerate
        # then use tail and cut to list only grnas
