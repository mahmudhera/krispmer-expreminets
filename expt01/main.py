import dna_jellyfish as jellyfish
import subprocess
import pandas as pd
import sys
from get_cfd_score import get_score
import trie
from itertools import chain, combinations, product
import sys
from multiprocessing import Process, Manager

num_threads = 48
num_cores = num_threads

def test_parallel(arg_key, arg_val, dict):
    dict[arg_key] = arg_val

def test():
    manager = Manager()
    processes = []
    dictionaries = []
    for i in range(num_cores):
        dict = manager.dict()
        key = i
        val = "val is " + str(i)
        dictionaries.append(dict)
        p = Process( target=test_parallel, args=(key, val, dict) )
        processes.append(p)
        p.start()

    for p in processes:
        p.join()
    for d in dictionaries:
        print(d)

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

def generate_inverted_specificity_from_genome_parallel(guides, qf_genome, qf_target, max_hd, target_count, dict):
    dict = {}
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
            dict[candidate] = 1.0 * val1 / (val2 * target_count)
        except:
            dict[candidate] = -1
            print(cutting_probability, get_score(reverse_complement(candidate), reverse_complement(mer)))

def generate_inverted_specificity_from_genome(guides, qf_genome, qf_target, max_hd = 3, target_count = 1):
    process_list = []
    dictionaries = []
    guides_per_thread = len(guides) / num_threads
    handled_so_far = 0
    manager = Manager()
    for i in range(num_cores):
        low_index = handled_so_far
        high_index = min( int( (i+1)*guides_per_thread ) , len(candidates) )
        guides_this_thread = guides[low_index:high_index]
        dict = manager.dict()
        p = Process(target=generate_inverted_specificity_from_genome_parallel, args = (guides, qf_genome, qf_target, max_hd, target_count, dict))
        proces_list.append(p)
        dictionaries.append(dict)
        p.start()

    for p in process_list:
        p.join()

    dict = {}
    for d in dictionaries:
        dict.update(d)
    return dict


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("error")
		print("usage: python main.py genome_jf_filename target_filename scores_filename max_hd out_filename <tgt_count, optional>")
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
