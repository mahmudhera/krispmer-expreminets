from matplotlib import pyplot as plt
import numpy as np
from collections import Counter
import pandas as pd
import sys

jitter = 0.06
x_label_values = [1, 2]

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("error")
        print("usage: python plot_boxes.py out_filename target_name genome_name")
        sys.exit(-1)

    filename = sys.argv[1]
    target_name = sys.argv[2]
    genome_name = sys.argv[3]

    scores_krispmer = pd.read_csv(filename, header=None, delimiter=' ')[2].to_list()
    scores_using_genome = pd.read_csv(filename, header=None, delimiter=' ')[1].to_list()

    scores = [scores_krispmer, scores_using_genome]
    labels = ['kRISP-meR', 'Genome']

    for i in range(2):
        y = scores[i]
        x = np.random.normal(x_label_values[i], jitter, len(y))
        plt.plot(x, y, 'x', alpha=0.2)
    plt.boxplot(scores, showfliers=False, widths=0.5)
    plt.ylabel('Inverse specificity')
    plt.xticks(x_label_values, labels, rotation=20)
    plt.title('Genome: ' + genome_name + ', target: ' + target_name)
    plt.savefig(filename + '.pdf')
