import dna_jellyfish as jellyfish
import subprocess
import pandas as pd
import sys
from get_cfd_score import get_score


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


def generate_inverted_specificity_from_genome(guides, qf_genome, qf_target, max_hd = 3, target_count = 1):
    """
    generates the inverted specificity score used in krispmer, but using the genome, not expectations
    :param guides: list of guides, which are 2-tuple of grna(23-nts) and strand('+' or '-')
    :param genome_jf_filename: jellyfish filename (generated for 23-mers)
    :param target_region_filename: the target string as fasta
    :return: dictionary containing the scores. dic[gRNA sequence] -> the inv_spec score
    """
    dic = {}
    for candidate in guides:
        trie = generate_adjacent_mers(candidate, max_hd)
        val1 = 0.0
        val2 = 0.0
        for mer in trie.keys():
            merDNA = jellyfish.MerDNA(mer)
            revcompMerDNA = jellyfish.MerDNA(reverse_complement(mer))
            cutting_probability = get_score(candidate, mer)
            val1 += max(qf[merDNA], qf[revcompMerDNA]) * cutting_probability
            val2 += max(tgt_qf[merDNA], tgt_qf[revcompMerDNA]) * cutting_probability
            #print (mer + ' ' + str(cutting_probability))
        try:
            dic[candidate] = 1.0 * val1 / (val2 * target_count)
        except:
            dic[candidate] = -1
    return dic


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("error")
		print("usage: python main.py genome_jf_filename target_filename scores_filename max_hd out_filename")
		sys.exit(-1)

	genome_jf_filename = sys.argv[1]
	target_filename = sys.argv[2]
	scores_filename = sys.argv[3]
	max_hd = int(sys.argv[4])
	out_filename = sys.argv[5]

	try:
		target_count = sys.argv[6]
	except:
		target_count = 1

	qf_target = generate_jf_file(target_filename)
	qf_genome = jellyfish.QueryMerFile(genome_jf_filename)

	grnas_in_positive = get_list_of_grna_strings(scores_filename)
	genome_scores = generate_inverted_specificity_from_genome(grnas_in_positive, qf_genome,
														qf_target, max_hd, target_count)
	print(genome_scores)
