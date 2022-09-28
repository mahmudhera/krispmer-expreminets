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
        genome_jf_fname = 'genome_counts.jf'
        scores_fname = 'scores_' + tgt_fname
        max_hd = '3'
        out_fname = 'tmt_scores'

        cmd_args = ['python', 'find_score_with_genome.py', genome_jf_fname, tgt_fname, scores_fname, max_hd, out_fname]
        subprocess.call(cmd_args)

        f = open(out_fname, 'r')
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            grna = line.split(' ')[0]
            score_with_genome = float(line.split(' ')[1])
            score_without_genome = float(line.split(' ')[2])
            list_grna_scores.append((score_with_genome, score_without_genome))

        break
    sys.stdout = open('all_scores', 'w')
    for score in list_grna_scores:
        print(score)
    sys.stdout.close()

if __name__ == "__main__":
    main()
