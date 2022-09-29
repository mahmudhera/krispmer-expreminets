import os
import subprocess
import sys

def main():
    list_grna_scores = []

    all_fasta_files = []
    for file in os.listdir('.'):
        if file.endswith('.fasta') and not file.startswith('scores') and not file.startswith('genome') and not file.startswith('target_NC'):
            all_fasta_files.append(file)

    print('Handling ' + str(len(all_fasta_files)) + ' targets...')
    c = 0
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
        c += 1
        print('Done with ' + str(c) + ' targtes.')

    sys.stdout = open('all_scores', 'w')
    for score in list_grna_scores:
        print(str(score[0]) + ' ' + str(score[1]))
    sys.stdout.close()

if __name__ == "__main__":
    main()
