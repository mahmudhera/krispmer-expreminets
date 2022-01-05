import dna_jellyfish as jellyfish
import subprocess
import pandas as pd

def get_off_target_count(target_jf_file, genome_jf_file, kmer_str):
	mer = jellyfish.MerDNA(kmer_str)
	cg1, ct1 = genome_jf_file[mer], target_jf_file[mer]
	mer.canonicalize()
	cg2, ct2 = genome_jf_file[mer], target_jf_file[mer]
	print(cg1 + cg2, ct1 + ct2)
	return max(0, cg1 + cg2 - ct1 - ct2)
	
def generate_jf_file(fasta_filename, jf_filename="temp"):
	# given a fasta file, use jellyfish to count 23-mers
	# returns a QueryMerFile
	jf_command = "jellyfish count -m 23 -s 100M -o " + jf_filename + " -t 8 -C " + fasta_filename
	args = jf_command.split(' ')
	subprocess.call(args)
	return jellyfish.QueryMerFile(jf_filename)
	
def get_list_of_grnas(scores_filename):
	l = pd.read_csv(scores_filename)[tgt_in_plus].to_list()
	print (l)

if __name__ == "__main__":
	genome_jf_filename = 'staphylococcus_genome.jf'
	test_kmer = 'CCAATTGGGGCCGTCTCTATAAT'
	target_qf = generate_jf_file("../../../data/staphylococcusAureus/target1")
	get_list_of_grnas("../../../data/staphylococcusAureus/scores1")
	
	genome_qf = jellyfish.QueryMerFile(genome_jf_filename)
	
	x = get_off_target_count(target_qf, genome_qf, test_kmer)
	print (x)
	
	print("compiles ok")