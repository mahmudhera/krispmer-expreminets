import dna_jellyfish as jellyfish

def get_off_target_count(target_jf_file, genome_jf_file, kmer):
	return 0

if __name__ == "__main__":
	target_jf_filename = 'test'
	genome_jf_filename = 'test'
	test_kmer = 'AGTCGTCGTCGTACGTGCGGG'
	
	target_qf = jellyfish.QueryMerFile(target_jf_filename)
	genome_qf = jellyfish.QueryMerFile(target_jf_filename)
	
	x = get_off_target_count(target_qf, genome_qf, test_kmer)
	print (x)
	
	print("compiles ok")