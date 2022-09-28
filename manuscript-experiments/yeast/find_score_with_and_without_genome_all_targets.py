import os
import subprocess

def main():
    list_grna_scores = []

    all_fasta_files = []
    for file in os.listdir('.'):
        if file.endswith('.fasta') and not file.startswith('scores'):
            all_fasta_files.append(file)

    for filename in all_fasta_files:
        tgt_fname = filename
        genome_jf_fname = 'mer_counts.jf'
        scores_fname = 'scores_' + tgt_fname
        max_hd = '3'
        out_fname = 'tmt_scores'

        cmd_args = ['python', 'find_score_with_genome.py', genome_jf_fname, tgt_fname, scores_fname, max_hd, out_fname]
        subprocess.call(cmd_args)

        print(tgt_fname, scores_fname)
        break

if __name__ == "__main__":
    main()
