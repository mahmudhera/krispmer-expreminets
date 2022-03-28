This is the first experiment.

Next meeting:
1. find the bug, fix it
1. count CHM13 kmers in reference
1. dig: which version of hg38 do they report using??

# Score with and without using ref. genome

The steps are as follows:

1. Run for multiple genomes for a hd=3
1. Generate jf file for the genome
1. Run commands like follows:
python main.py ../runs-on-genomes/staphylococcus-aureus/mer_counts.jf ../runs-on-genomes/staphylococcus-aureus/target4.fasta ../runs-on-genomes/staphylococcus-aureus/scores4 3 staphylococcus-aureus/out_for_target4


To do the plotting:
conda activate plotting
