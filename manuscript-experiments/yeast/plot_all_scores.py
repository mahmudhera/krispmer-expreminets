from matplotlib import pyplot as plt
import pandas as pd

scores_with_genomes = []
scores_without_genomes = []

f = open('all_scores', 'r')
lines = f.readlines()
f.close()
for line in lines:
    scores_with_genomes.append( float(line.split(' ')[0]) )
    scores_without_genomes.append( float(line.split(' ')[1]) )

plt.scatter(scores_with_genomes, scores_without_genomes, s=5)
plt.plot( [min(scores_with_genomes), max(scores_with_genomes)], [min(scores_with_genomes), max(scores_with_genomes)], linestyle='--' )
#plt.xlim(0,2)
#plt.ylim(0,5)
plt.xlabel('Score with full genome')
plt.ylabel('krispmer score without genome')
plt.savefig('scores_for_all_targets.pdf')
