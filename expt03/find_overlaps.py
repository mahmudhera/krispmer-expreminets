import pandas as pd
from sys import argv
import dna_jellyfish as jellyfish
import subprocess

cut_off_score = 1.5

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

def get_guides_from_kr_file(fname):
    df = pd.read_csv(fname, delimiter=',')
    positive_grnas = df[ df['inverse_specificity']<=cut_off_score ]['tgt_in_plus'].tolist()
    negative_grnas = df[ df['inverse_specificity']<=cut_off_score ]['tgt_in_minus'].tolist()
    return positive_grnas, negative_grnas

def get_guides_from_gs_file(fname):
    df = pd.read_excel(fname)
    col_name = 'Number of off-targets'
    df = df[ df[col_name] == 0 ]
    return df['Target-Seq'].tolist()

def find_common(gs_grnas, kr_grnas_in_pos, kr_grnas_in_neg):
    count = 0
    for guide in gs_grnas:
        found = False
        for guide_pos in kr_grnas_in_pos:
            if guide in guide_pos:
                count+=1
                found = True
                break
        if found:
            continue
        for guide_neg in kr_grnas_in_neg:
            if guide in guide_neg:
                count+=1
                found=True
                break
    return count

def generate_jf_file(fasta_filename, jf_filename="temp"):
	jf_command = "jellyfish count -m 23 -s 100M -o " + jf_filename + " -t 8 -C " + fasta_filename
	args = jf_command.split(' ')
	subprocess.call(args)
	return jellyfish.QueryMerFile(jf_filename)

def get_off_target_count(target_jf_file, genome_jf_file, kmer_str):
	mer = jellyfish.MerDNA(kmer_str)
	cg1, ct1 = genome_jf_file[mer], target_jf_file[mer]
	mer = jellyfish.MerDNA( reverse_complement(kmer_str) )
	cg2, ct2 = genome_jf_file[mer], target_jf_file[mer]
	print (cg1+cg2, ct1+ct2)
	return max(0, cg1 + cg2 - ct1 - ct2)

def get_off_target_counts_krispmer(grnas_in_positive):
    grnas = grnas_in_positive
    off_tgt_counts = []
    for grna in grnas:
        off_tgt_counts.append( get_off_target_count(grna[:-3]) )
    return off_tgt_counts

def get_off_target_counts_guidescan(grnas):
    off_tgt_counts = []
    for grna in grnas:
        off_tgt_counts.append( get_off_target_count(grna) )
    return off_tgt_counts

if __name__=="__main__":
    print("Usage: python find_overlaps.py krispmer_fname guidescan_fname genome_jf_fname target_fname")
    kr_fname = argv[1]
    gs_fname = argv[2]
    genome_jf_fname = argv[3]
    generate_jf_file(argv[4])
    target_jf_fname = "temp"

    grnas_in_pos, grnas_in_neg = get_guides_from_kr_file(kr_fname)
    print("Number of grnas by krispmer:")
    print(len(grnas_in_neg))
    gs_grnas = get_guides_from_gs_file(gs_fname)
    print("Number of grnas by guidescan:")
    print(len(gs_grnas))
    common = find_common(gs_grnas, grnas_in_pos, grnas_in_neg)
    print("Only by krispmer \tcommon\tOnly by guidescan")
    print(str(len(grnas_in_neg) - common) + '\t\t\t' + str(common) + '\t' + str(len(gs_grnas) - common))
    print('')

    print( get_off_target_counts_krispmer(grnas_in_pos) )
    print( get_off_target_counts_guidescan(gs_grnas) )
