import dna_jellyfish as jellyfish

def get_off_target_count(target_jf_file, genome_jf_file, kmer_str):
	mer = jellyfish.MerDNA(kmer_str)
	mer = mer.canonicalize()
	return max(0, genome_jf_file[mer]-target_jf_file[mer])

if __name__ == "__main__":
	target_jf_filename = 'mer_counts.jf'
	genome_jf_filename = 'staphylococcus_genome.jf'
	test_kmer = 'AGTCGTCGTCGTACGTGCGGG'
	
	target_qf = jellyfish.QueryMerFile(target_jf_filename)
	genome_qf = jellyfish.QueryMerFile(target_jf_filename)
	
	x = get_off_target_count(target_qf, genome_qf, test_kmer)
	print (x)
	
	print("compiles ok")