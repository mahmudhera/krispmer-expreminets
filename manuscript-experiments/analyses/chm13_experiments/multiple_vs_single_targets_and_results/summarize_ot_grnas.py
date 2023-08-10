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
#methods = ['crispor', 'guidescan2', 'krispmer', 'krispmer2.0', 'krispmer3.0']
methods = ['crispor', 'guidescan2', 'krispmer']
num_targets = 3
ot_categories = ['0 mismatch OT', '1 mismatch OT', '2/more mismatch OT']
transcript_ids = ['ENST00000375050.6', 'ENST00000568530.5', 'ENST00000261405.10']

for target_id in range(num_targets):

    off_target_counts_0_mm = list([0]*len(methods))
    off_target_counts_1_mm = list([0]*len(methods))
    off_target_counts_2_or_more_mm = list([0]*len(methods))

    methods_to_num_occurrences_dict = {}

    for i in range( len(methods) ):
        method_name = methods[i]
        aggregated_filename = f'{file_str}{method_name}'
        df = pd.read_csv(aggregated_filename, delimiter=' ', header=None, names=['target_name', 'grna', 'a', 'b', 'c', 'num_occurrences_of_grna_in_genome'])
        df = df[ df['target_name']==f'target_{target_id+1}.fasta' ]
        targets = df['target_name'].tolist()
        grnas = df['grna'].tolist()
        a_s = df['a'].tolist()
        b_s = df['b'].tolist()
        c_s = df['c'].tolist()
        num_occurrences_of_grna_in_genome_list = df['num_occurrences_of_grna_in_genome'].tolist()
        methods_to_num_occurrences_dict[method_name] = num_occurrences_of_grna_in_genome_list

        for target_name, grna, a, b, c, num_occurrences_of_grna_in_genome in list(zip(targets, grnas, a_s, b_s, c_s, num_occurrences_of_grna_in_genome_list)):
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
    plt.title(f'gRNAs identified by tools for transcript id {transcript_ids[target_id]}')
    plt.savefig(pdf_filename)
    plt.clf()

    pdf_filename = f'num_occurrences_plot_target_{target_id+1}.pdf'
    fig, ax = plt.subplots()
    ax.boxplot(methods_to_num_occurrences_dict.values())
    ax.set_xticklabels(methods_to_num_occurrences_dict.keys())
    plt.xlabel('gRNA finding tools')
    plt.ylabel('Number of occurrences of the gRNAs in CHM13')
    plt.title(f'Number of times the gRNAs appear in CHM13 for transcript: {transcript_ids[target_id]}')
    plt.savefig(pdf_filename)
    plt.clf()
