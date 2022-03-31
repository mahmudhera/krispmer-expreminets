from matplotlib import pyplot as plt
import numpy as np
from collections import Counter
import pandas as pd
import sys
from scipy.signal import savgol_filter
from scipy.stats.stats import pearsonr

jitter = 0.06
x_label_values = [1, 2]

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("error")
        print("usage: python <py filename> out_filename target_name genome_name")
        sys.exit(-1)

    filename = sys.argv[1]
    target_name = sys.argv[2]
    genome_name = sys.argv[3]

    scores_krispmer = pd.read_csv(filename, header=None, delimiter=' ')[2].to_list()
    scores_using_genome = pd.read_csv(filename, header=None, delimiter=' ')[1].to_list()

    plt.scatter( scores_krispmer, scores_using_genome )

    plt.suptitle('Genome: ' + genome_name + ', target: ' + target_name)
    plt.title('Correlation coefficient: %0.4lf' % pearsonr(scores_krispmer, scores_using_genome)[0], fontsize=8)
    plt.savefig(filename + '_x_y.pdf')
