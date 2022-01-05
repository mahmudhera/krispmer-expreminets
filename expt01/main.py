import dna_jellyfish as jellyfish

def get_off_target_count(target_jf_file, genome_jf_file, kmer_str):
	mer = jellyfish.MerDNA(kmer_str)
	print(genome_jf_file[mer], target_jf_file[mer])
	mer.canonicalize()
	print(genome_jf_file[mer], target_jf_file[mer])
	return max(0, genome_jf_file[mer]-target_jf_file[mer])
	
#def generate_jf_file(fasta_filename):
	# given a fasta file, use jellyfish to count 23-mers
	# returns a QueryMerFile

if __name__ == "__main__":
	target_jf_filename = 'mer_counts.jf'
	genome_jf_filename = 'staphylococcus_genome.jf'
	test_kmer = 'CCAATTGGGGCCGTCTCTATAAT'
	
	target_qf = jellyfish.QueryMerFile(target_jf_filename)
	genome_qf = jellyfish.QueryMerFile(target_jf_filename)
	
	x = get_off_target_count(target_qf, genome_qf, test_kmer)
	print (x)
	
	print("compiles ok")