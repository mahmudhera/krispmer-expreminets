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

    off_targets = [score-1 for score in scores_using_genome]

    scores_ot = list(zip(scores_krispmer, off_targets))
    scores_ot.sort(key=lambda x: x[0])

    cdf = []
    cumulative_off_target = 0
    for (score, off_target) in scores_ot:
        cumulative_off_target += off_target
        cdf.append( (score, cumulative_off_target) )
        if score > 2.0:
            break

    xs = [tuple[0] for tuple in cdf]
    ys = [tuple[1] for tuple in cdf]
    ys = savgol_filter(ys, 3, 2)

    plt.plot( xs, ys )

    plt.title('Genome: ' + genome_name + ', target: ' + target_name)
    plt.savefig(filename + '.pdf')
