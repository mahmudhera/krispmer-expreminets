# THE script that generates plots

import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np
import os

# Step 2: Choose the font 'Times New Roman'
title_font = {'fontname': 'Arial', 'fontsize': 12}
axis_label_font = {'fontname': 'Arial', 'fontsize': 10}
tick_label_font = {'fontname': 'Arial', 'fontsize': 10}

# Step 3: Set the font properties using rcParams
plt.rcParams['font.family'] = 'sans-serif'  # Set the default font family
plt.rcParams['font.sans-serif'] = 'Arial'  # Set the default font to Times New Roman

file_str = 'ot_summary_chm13_'
#methods = ['crispor', 'guidescan2', 'krispmer', 'krispmer2.0', 'krispmer3.0']
methods = ['crispor', 'guidescan2', 'krispmer']
ot_categories = ['0 mismatch OT', '1 mismatch OT', '2/more mismatch OT']

for target_filename in os.listdir('targets'):
    off_target_counts_0_mm = list([0]*len(methods))
    off_target_counts_1_mm = list([0]*len(methods))
    off_target_counts_2_or_more_mm = list([0]*len(methods))

    for i in range( len(methods) ):
        method_name = methods[i]
        aggregated_filename = f'{file_str}{method_name}'
        df = pd.read_csv(aggregated_filename, delimiter=' ', header=None, names=['target_name', 'grna', 'a', 'b', 'c'])
        df = df[ df['target_name']==target_filename ]
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

    splitted = target_filename.split('.')[0]
    pdf_filename = f'ot_summary_pdfs/ot_summary_{splitted}.pdf'
    fig, ax = plt.subplots()
    ax.bar(methods, off_target_counts_0_mm, label='0 mismatch to OTs')
    ax.bar(methods, off_target_counts_1_mm, bottom=off_target_counts_0_mm, label='1 mismatch to OTs')
    ax.bar(methods, off_target_counts_2_or_more_mm, bottom=np.add(off_target_counts_0_mm, off_target_counts_1_mm), label='2/more mismatches to OTs')
    ax.legend(loc='upper right')
    plt.xlabel('gRNA finding tools')
    plt.ylabel('Number of gRNAs identified')
    splitted = target_filename.split('_')[-1].split('.')[0]
    plt.title(f'gRNAs identified by tools for target: {splitted}')
    plt.savefig(pdf_filename)
    plt.clf()
