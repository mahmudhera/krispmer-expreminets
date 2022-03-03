import dna_jellyfish as jellyfish
import subprocess
import pandas as pd
import sys

# number of off-target counts of kmer_str in genome
# requires two jf files
def get_off_target_count(target_jf_file, genome_jf_file, kmer_str):
	mer = jellyfish.MerDNA(kmer_str)
	cg1, ct1 = genome_jf_file[mer], target_jf_file[mer]
	mer.canonicalize()
	cg2, ct2 = genome_jf_file[mer], target_jf_file[mer]
	print (cg1+cg2, ct1+ct2)
	return max(0, cg1 + cg2 - ct1 - ct2)

# given a fasta file, use jellyfish to count 23-mers
# returns a QueryMerFile
def generate_jf_file(fasta_filename, jf_filename="temp"):
	jf_command = "jellyfish count -m 23 -s 100M -o " + jf_filename + " -t 8 -C " + fasta_filename
	args = jf_command.split(' ')
	subprocess.call(args)
	return jellyfish.QueryMerFile(jf_filename)

# returns grnas in +ve strand
def get_list_of_grna_strings(scores_filename):
	return pd.read_csv(scores_filename)['tgt_in_plus'].to_list()



if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("error")
		print("usage: python main.py working_direcory genome_jf_filename target_filename scores_filename max_hd")
		sys.exit(-1)

	folder_name = sys.argv[1]
	genome_jf_filename = sys.argv[2]
	target_filename = sys.argv[3]
	scores_filename = sys.argv[4]
	max_hd = int(sys.argv[5])

	print(folder_name, target_filename)

	#qf_target = generate_jf_file()
	# create two jf files
	# copy fn from github
	# call that
