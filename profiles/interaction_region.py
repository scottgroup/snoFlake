#!/usr/bin/env python

""" interaction_region.py: plot the interaction region profile in the snoRNAs """


import pandas as pd
import sys
import os
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO
import re
import numpy as np
import seaborn as sns
import subprocess

sns.set_context('talk')


def identify_boxes_CD(sno_seq):
    """
    Identify the position of the boxes C, D' and D
    :param sno_seq: sequence of the snoRNA
    :return: position of the boxes C, D' and D
    """
    sno_seq = sno_seq.upper()
    group_boxes = re.search('^([GTCNA]+?)?([ACGT][GAT]GA[TG]G[ATG])[AGTNC]+$', sno_seq)
    try:
        ext5 = group_boxes.group(1)
        posC = len(ext5)
    except(TypeError, AttributeError):
        posC = None

    group_boxes = re.search('^[AGTNC]+(CTGA)([GTCNA]+?)?$', sno_seq)
    try:
        ext3 = group_boxes.group(2)
        len_ext3 = len(ext3)
        posD = len(sno_seq) - len_ext3 - 4
    except(TypeError, AttributeError):
        posD = None

    if posC is None:
        substart = 0
    else:
        #if box C is identified, only look for box D' sequence after the box C
        substart = posC + 7
    if posD is None:
        subend = None
    else:
        # if box D is identified, only look for box D' sequence before the box D
        subend = posD
    delta = substart
    subgroup_boxes = re.search('^([AGTNC]+?)?(CTGA)[GTCNA]+?$', sno_seq[substart: subend])
    try:
        pos_Dp = len(subgroup_boxes.group(1))
        pos_Dp += delta
    except(TypeError, AttributeError):
        pos_Dp = None

    return posC, pos_Dp, posD


def plot_region(df, snoid, sno_entry, score, nb_windows, fpath, relplot_path):
    df = df.groupby(['sno_window_start', 'sno_window_end'])['snoid'].count().reset_index()
    sno_seq = sno_entry['seq']
    # compute snoRNA folding using RNAplot and the entropy of the folding using relplot
    temp_fasta = os.path.join(fpath, 'sno_int_region', snoid + '.fa')
    with open(temp_fasta, 'w') as w:
        w.write('>' + snoid + '\n' + str(sno_seq))
    cmd_rf = ['RNAfold', '-p', temp_fasta]
    res2 = subprocess.check_output(cmd_rf)
    pairing = res2.decode().split()[2]
    rel_file = os.path.join(fpath, 'sno_int_region', snoid + '_rel.ps')
    with open(rel_file, 'w') as rel_w:
        cmd_relplot = [relplot_path,  os.path.join(fpath, 'sno_int_region', snoid + '_ss.ps'),
                       os.path.join(fpath, 'sno_int_region', snoid + '_dp.ps')]
        res3 = subprocess.run(cmd_relplot, stdout=rel_w)
    entropy = []
    line = ''
    with open(rel_file) as rel_r:
        while line.startswith('/S') is False:
            line = rel_r.readline()
        for line in rel_r:
            line = line.strip()
            if line == '] def':
                break
            entropy.append(float(line))
    os.remove(temp_fasta)

    # identify the position of the boxes
    pos_boxes = identify_boxes_CD(str(sno_seq))

    # count the number of interaction for each nucleotide
    pos = range(len(sno_seq))
    zero = [0] * len(pos)
    init_count = [[snoid, pos[i], zero[i]] for i in range(len(pos))]
    df_init = pd.DataFrame(init_count, columns=['snoid', 'pos', 'count'])
    for idx, row in df.iterrows():
        df_init.loc[(df_init.pos >= row['sno_window_start']) &
                    (df_init.pos < row['sno_window_end']), 'count'] += row['snoid']

    # plot the number of interaction per nucleotide, the folding and the entropy
    fig, axes = plt.subplots(2, 2, sharex='col', figsize=[15,10],
                             gridspec_kw={'width_ratios': [100, 1], 'height_ratios': [15, 1]})
    ax = axes[0,0]
    ax2 = axes[1,0]
    axes[0, 1].remove()
    ax.set_title('%s interaction regions with at least %d consecutive windows having a probability >= %s%%' % (
        snoid, nb_windows, score))
    ax.set_ylabel('Nb of predicted interactions')
    ax.step(x=df_init.pos, y=df_init['count'], where='post')
    sno_seq = str(sno_seq).upper()
    sno_label = [n + '\n' + pairing[i] for i, n in enumerate(sno_seq)]
    y_lim = ax.get_ylim()
    boxc_color = 'mediumseagreen'
    boxd_color = 'mediumorchid'
    for b, p in enumerate(pos_boxes):
        if p is None:
            continue
        if b == 0:
            l = 7
            color = boxc_color
        else:
            l = 4
            color = boxd_color
        ax.fill_between(range(p, p + l), -10000, df_init['count'].max() + 10000, facecolor=color, alpha=0.2)
    legend_elements = [Patch(facecolor=boxc_color, label='Box C', alpha=0.2),
                       Patch(facecolor=boxd_color, label='Box D\'/D', alpha=0.2)]

    im = ax2.imshow([entropy], vmin=min(entropy), vmax=max(entropy), cmap='gist_rainbow', aspect='auto')
    ax3 = axes[1,1]
    plt.colorbar(im, cax=ax3, ax=ax2)
    ax3.set_ylabel('Positional\nentropy')

    ax.set_ylim(y_lim)
    ax.set_xticks(pos)
    ax.set_xticklabels(sno_label)
    ax.xaxis.set_tick_params(which='both', labelbottom=True)
    plt.xticks(fontsize=12)
    ax.legend(handles=legend_elements)
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.tight_layout()
    for spine in ax2.spines.values():
        spine.set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)

    fig.savefig(os.path.join(fpath, 'sno_int_region', '%s_interaction_region.%s_%d.svg' % (snoid, score, nb_windows)))
    plt.close(fig)


def plot_all_sno(df_all, score, nb_windows, fpath, ytype):
    pos = np.arange(0.0000, 1.0001, 0.0001)
    zero = [0] * len(pos)
    init_count = [['all', pos[i], zero[i]] for i in range(len(pos))]
    df_init = pd.DataFrame(init_count, columns=['proportion', 'pos', 'count'])
    for idx, row in df_all.iterrows():
        df_init.loc[(df_init.pos >= row['sno_window_start']) &
                    (df_init.pos < row['sno_window_end']), 'count'] += row['proportion']
    fig = plt.figure(figsize=[15, 8])
    ax = fig.add_subplot()
    ax.set_title('All snoRNA interaction regions with at least %d consecutive windows having a probability >= %s%%' % (
        nb_windows, score))
    if ytype == 'count':
        ax.set_ylabel('Nb of predicted interactions')
    elif ytype == 'prop':
        ax.set_ylabel('Normalized nb of predicted interactions for each snoRNA')
    ax.set_xlabel('Relative position')
    ax.step(x=df_init.pos, y=df_init['count'], where='post')
    fig.tight_layout()
    fig.savefig(os.path.join(fpath, 'sno_int_region', 'ALL_interaction_region.%s.%s_%d.svg' % (ytype, score, nb_windows)))
    plt.close(fig)


def main():
    inpath = os.path.abspath(sys.argv[1]) # path to the interaction files
    sno_fasta = os.path.abspath(sys.argv[2]) # fasta file of snoRNA sequences
    score = sys.argv[3] # score cutoff (used for filename and fig title)
    nb_windows = int(sys.argv[4]) # nb of consecutive windows cutoff
    relplot_path = os.path.abspath(sys.argv[5]) # path to relplot utility from viennaRNA

    # create folder for figures
    os.makedirs(os.path.join(inpath, 'sno_int_region'), exist_ok=True)
    os.chdir(os.path.join(inpath, 'sno_int_region'))

    # read snoRNA fasta
    sno_dict = SeqIO.to_dict(SeqIO.parse(sno_fasta, 'fasta'))
    sno_dict = {k.split('|')[0]: {'chromo': k.split('|')[1], 'seq': v.seq} for k, v in sno_dict.items()}

    df_all = pd.DataFrame()
    for snoid in sno_dict.keys():
        df_sno = pd.DataFrame()
        if snoid in ['NR_003072', 'NR_003073', 'ENSG00000252284']:
            continue
        infile = os.path.join(inpath, 'pred_%s.%s_%s.gene.tsv' % (snoid, score, nb_windows))
        try:
            df = pd.read_csv(infile,
                             sep='\t',
                             names=['seqname', 'target_window_start', 'target_window_end',
                                    'snoid', 'count', 'strand', 'sno_window_start',
                                    'sno_window_end', 'mean_score', 'min_score', 'max_score', 'target'])
            df = df[df.snoid != 'snoid']  # drop header if present
            df[['count', 'sno_window_start', 'sno_window_end']] = df[
                ['count', 'sno_window_start', 'sno_window_end']].astype(float)
            df = df.drop(['target'], axis=1)
            df = df.drop_duplicates()
            df = df[['snoid', 'sno_window_start', 'sno_window_end', 'count']]
            df = df[df['count'] >= nb_windows]
            df_sno = df_sno.append(df)
        except FileNotFoundError:
            print('File not found: %s' % infile, file=sys.stderr)
            continue
        if df_sno.empty:
            print('No %s interaction met the criteria: %s, %d'%(snoid, score, nb_windows), file=sys.stderr)
            continue
        plot_region(df_sno, snoid, sno_dict[snoid], score, nb_windows, inpath, relplot_path)
        df_sno['length'] = len(sno_dict[snoid]['seq'])
        df_sno['proportion'] =  1 / len(df_sno)
        df_all = df_all.append(df_sno)
    df_all['sno_window_start'] = ((df_all['sno_window_start']) / df_all['length']).round(4)
    df_all['sno_window_end'] = ((df_all['sno_window_end']) / df_all['length']).round(4)
    df_count = df_all.groupby(['sno_window_start', 'sno_window_end'])['proportion'].count().reset_index()
    df_prop = df_all.groupby(['sno_window_start', 'sno_window_end'])['proportion'].sum().reset_index()
    # with "prop" option, the profiles are normalize by the number of interactions for each snoRNA, the total of all
    # interactions for one snoRNA is 1.
    plot_all_sno(df_prop, score, nb_windows, inpath, 'prop')
    # with "count" option, all the interactions are counted once.
    plot_all_sno(df_count, score, nb_windows, inpath, 'count')


if __name__ == '__main__':
    main()
