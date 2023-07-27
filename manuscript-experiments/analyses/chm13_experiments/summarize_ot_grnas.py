import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np

# Step 2: Choose the font 'Times New Roman'
title_font = {'fontname': 'Arial', 'fontsize': 12}
axis_label_font = {'fontname': 'Arial', 'fontsize': 10}
tick_label_font = {'fontname': 'Arial', 'fontsize': 10}

# Step 3: Set the font properties using rcParams
plt.rcParams['font.family'] = 'sans-serif'  # Set the default font family
plt.rcParams['font.sans-serif'] = 'Arial'  # Set the default font to Times New Roman

file_str = 'summary_'
methods = ['crispor', 'guidescan2', 'krispmer']
num_targets = 3
ot_categories = ['0 mismatch OT', '1 mismatch OT', '2/more mismatch OT']


for target_id in range(num_targets):

    off_target_counts_0_mm = [0, 0, 0]
    off_target_counts_1_mm = [0, 0, 0]
    off_target_counts_2_or_more_mm = [0, 0, 0]

    for i in range( len(methods) ):
        method_name = methods[i]
        aggregated_filename = f'{file_str}{method_name}'
        df = pd.read_csv(aggregated_filename, delimiter=' ', header=None, names=['target_name', 'grna', 'a', 'b', 'c'])
        df = df[ df['target_name']==f'target_{target_id+1}.fasta' ]
        targets = df['target_name'].tolist()
        grnas = df['grna'].tolist()
        a_s = df['a'].tolist()
        b_s = df['b'].tolist()
        c_s = df['c'].tolist()

        for target_name, grna, a, b, c in list(zip(targets, grnas, a_s, b_s, c_s)):
            if a+b+c == 0:
                off_target_counts_2_or_more_mm[i] += 1
                continue
            if a > 0:
                off_target_counts_0_mm[i] += 1
                continue
            if b > 0:
                off_target_counts_1_mm[i] += 1
                continue
            if c > 0:
                off_target_counts_2_or_more_mm[i] += 1

    pdf_filename = f'ot_summary_plot_target_{target_id+1}.pdf'
    fig, ax = plt.subplots()
    ax.bar(methods, off_target_counts_0_mm, label='0 mismatch to OTs')
    ax.bar(methods, off_target_counts_1_mm, bottom=off_target_counts_0_mm, label='1 mismatch to OTs')
    ax.bar(methods, off_target_counts_2_or_more_mm, bottom=np.add(off_target_counts_0_mm, off_target_counts_1_mm), label='2/more mismatches to OTs')
    ax.legend(loc='upper right')
    plt.xlabel('gRNA finding tools')
    plt.ylabel('Number of gRNAs identified')
    plt.savefig(pdf_filename)
    plt.clf()
