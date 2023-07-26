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

def plot_clustered_stacked(dfall, labels=None, title="multiple stacked bar plot",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot.
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""
    plt.figure(figsize=(5.5, 5.5))
    n_df = len(dfall)
    n_col = len(dfall[0].columns)
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    hatch_styles = ['', '//', '\\\\']
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.edgecolor='black'
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                #rect.set_hatch(H * int(i / n_col)) #edited part
                rect.set_hatch(hatch_styles[int(i / n_col)]) #edited part
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 30)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=hatch_styles[int(i)]))

    l1 = axe.legend(h[:n_col], l[:n_col], loc='upper right')
    if labels is not None:
        l2 = plt.legend(n, labels, loc='upper left')
    axe.add_artist(l1)
    return axe

organisms = ['human', 'mouse', 'yeast']
file_str = 'off_target_summary_'
methods = ['crispor', 'guidescan', 'krispmer']
types = ['protein', 'noncoding', 'duplicate']

# 1 to 50 - noncoding, 51-60 - protein, target_NC: duplicate masked

records = {}
# records[organism][type][method] = count

for organism in organisms:
    records[organism] = {}
    for type in types:
        records[organism][type] = {}
        for method in methods:
            records[organism][type][method] = [0,0,0,0]

for organism in organisms:
    for method in methods:
        aggregated_filename = f'{organism}/{file_str}{method}'
        df = pd.read_csv(aggregated_filename, delimiter=' ', header=None, names=['target_name', 'grna', 'a', 'b', 'c'])
        targets = df['target_name'].tolist()
        grnas = df['grna'].tolist()
        a_s = df['a'].tolist()
        b_s = df['b'].tolist()
        c_s = df['c'].tolist()

        type = None
        for target_name, grna, a, b, c in list(zip(targets, grnas, a_s, b_s, c_s)):
            if organism == 'yeast':
                if 'target_NC' in target_name:
                    type = 'duplicate'
                else:
                    num = int(target_name.split('_')[0][6:])
                    if num <= 50:
                        type = 'protein'
                    else:
                        type = 'noncoding'
            else:
                if 'gene' in target_name:
                    type = 'protein'
                if 'intergen' in target_name:
                    type = 'noncoding'
                if 'duplicate' in target_name:
                    type = 'duplicate'

            if a+b+c == 0:
                records[organism][type][method][0] += 1
                records[organism][type][method][3] += 1
                continue
            if a > 0:
                records[organism][type][method][1] += 1
                continue
            if b > 0:
                records[organism][type][method][2] += 1
                continue
            if c > 0:
                records[organism][type][method][3] += 1

ot_categories = ['0 mismatch OT', '1 mismatch OT', '2/more mismatch OT']

description = {
    'protein' : 'protein coding',
    'duplicate' : 'duplicate masked',
    'noncoding' : 'protein non-coding'
}

for type in types:
    # start a figure
    df_list = []

    for method in methods:
        df_data = []
        for organism in organisms:
            df_data.append( list(records[organism][type][method])[1:] )
        df = pd.DataFrame( df_data, index=organisms, columns=ot_categories )
        df_list.append(df)

    plot_clustered_stacked(df_list, methods, title=f'Number of gRNAs with Off-target counts in {description[type]} regions')
    plt.tight_layout()
    plt.savefig(f'{type}_ot_summary.pdf')
    plt.clf()
