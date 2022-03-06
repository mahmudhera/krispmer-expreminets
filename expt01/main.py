import dna_jellyfish as jellyfish
import subprocess
import pandas as pd
import sys
from get_cfd_score import get_score
import trie
from itertools import chain, combinations, product
import sys

def complement(seq):
    """
    generates the complement sequence, e.g.: ACGT-->TCGA
    :param seq: dna sequence
    :return: complement
    """
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(str(seq))
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)

def reverse_complement(s):
    return complement(s[::-1])

def hamming_circle(s, n, alphabet, trie):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.
    >>> sorted(hamming_circle('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_circle('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_circle('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']
    """
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            trie[''.join(cousin)]=1

def hamming_ball(s, n, alphabet, trie):
	"""Generate strings over alphabet whose Hamming distance from s is
	less than or equal to n.
	>>> sorted(hamming_ball('abc', 0, 'abc'))
	['abc']
	>>> sorted(hamming_ball('abc', 1, 'abc'))
	['aac', 'aba', 'abb', 'abc', 'acc', 'bbc', 'cbc']
	>>> sorted(hamming_ball('aaa', 2, 'ab'))
	['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba']
	"""
	for i in range(n+1):
		hamming_circle(s, i, alphabet, trie)

def generate_adjacent_mers(sequence, max_hamming_distance):
	alphabet = 'AGCT'
	t = trie.trie()
	hamming_ball(sequence, max_hamming_distance, alphabet, t)
	return t

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
	plus_ = pd.read_csv(scores_filename)['tgt_in_plus'].to_list()
	minus_ = pd.read_csv(scores_filename)['tgt_in_minus'].to_list()
	strand_ = pd.read_csv(scores_filename)['strand'].to_list()
	grnas = []
	for (p, m, s) in list(zip(plus_, minus_, strand_)):
		if s == '+':
			grnas.append(p)
		else:
			grnas.append(m)
	return grnas

# returns dictionary[grna]->our_score
def get_krispmer_scores(scores_filename):
	plus_ = pd.read_csv(scores_filename)['tgt_in_plus'].to_list()
	scores_ = pd.read_csv(scores_filename)['inverse_specificity'].to_list()
	d = {}
	for (p,s) in list(zip(plus_, scores_)):
		d[p] = s
	return d

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
            #print(cutting_probability, get_score(reverse_complement(candidate), reverse_complement(mer)))
            val1 += (qf_genome[merDNA] + qf_genome[revcompMerDNA]) * cutting_probability
            val2 += (qf_target[merDNA] + qf_target[revcompMerDNA]) * cutting_probability
            #print (mer + ' ' + str(cutting_probability))
        try:
            dic[candidate] = 1.0 * val1 / (val2 * target_count)
        except:
            dic[candidate] = -1
            print(cutting_probability, get_score(reverse_complement(candidate), reverse_complement(mer)))
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
	krispmer_scores = get_krispmer_scores(scores_filename)

	sys.stdout = open(out_filename, 'w')
	for grna in genome_scores.keys():
		try:
			print(grna, genome_scores[grna], krispmer_scores[grna])
		except:
			print(grna, genome_scores[grna], krispmer_scores[reverse_complement(grna)])
	sys.stdout.close()
