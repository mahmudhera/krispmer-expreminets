import dna_jellyfish as jellyfish
import subprocess
import pandas as pd
import sys

def get_off_target_count(target_jf_file, genome_jf_file, kmer_str):
	mer = jellyfish.MerDNA(kmer_str)
	cg1, ct1 = genome_jf_file[mer], target_jf_file[mer]
	mer.canonicalize()
	cg2, ct2 = genome_jf_file[mer], target_jf_file[mer]
	print (cg1+cg2, ct1+ct2)
	return max(0, cg1 + cg2 - ct1 - ct2)

def generate_jf_file(fasta_filename, jf_filename="temp"):
	# given a fasta file, use jellyfish to count 23-mers
	# returns a QueryMerFile
	jf_command = "jellyfish count -m 23 -s 100M -o " + jf_filename + " -t 8 -C " + fasta_filename
	args = jf_command.split(' ')
	subprocess.call(args)
	return jellyfish.QueryMerFile(jf_filename)

def get_list_of_grna_strings(scores_filename):
	return pd.read_csv(scores_filename)['tgt_in_plus'].to_list()

if __name__ == "__main__":
	if len(sys.argv) < 10:
		print("error")
		print("usage: python main.py working_direcory genome_jf_filename target_filename scores_filename max_hd")
		sys.exit(-1)

	print(sys.argv)

	folder_name = sys.argv[2]
	genome_jf_filename = sys.argv[3]
	target_filename = sys.argv[4]
	scores_filename = sys.argv[5]
	max_hd = int(sys.argv[6])
