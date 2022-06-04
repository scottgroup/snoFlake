#!/usr/bin/env python

""" Adpated from Gabrielle's script """
""" nts_flanking_exons.py: plot the interaction profiles inside the exons and in the 100 nts upstream and downstream """


import pandas as pd
import sys
import os
from pybedtools import BedTool
from matplotlib import pyplot as plt
import numpy as np
from multiprocessing import Pool
import seaborn as sns

sns.set_context('talk')


def normalize_count_up_down(df, len_dict, norm_factor):
    df = df[['rel_pos', 'counts']].groupby('rel_pos').counts.sum().reset_index()
    df['nb_introns'] = df.rel_pos.map(len_dict)
    df['nb_introns'] = df['nb_introns'].fillna(0)
    df['nb_introns'] = df['nb_introns'][::-1].cumsum()[::-1]
    df['counts'] = df['counts'] / df['nb_introns'] * norm_factor
    df = df.drop(['nb_introns'], axis=1)
    return df


def prep_up(df_up, len_dict, norm_factor):
    df_up.loc[df_up.strand == '+', 'rel_pos'] -= 101
    df_up['rel_pos'] = df_up['rel_pos'].abs()
    df_up = normalize_count_up_down(df_up, len_dict['up'], norm_factor)
    df_up['rel_pos'] = -df_up['rel_pos'] + 1
    return df_up


def prep_down(df_down, len_dict, norm_factor):
    df_down.loc[df_down.strand == '-', 'rel_pos'] = 101 - df_down['rel_pos']
    df_down = normalize_count_up_down(df_down, len_dict['down'], norm_factor)
    df_down['rel_pos'] -= 1
    return df_down


def prep_exon(df_exon, nb_exons, norm_factor):
    df_exon['length'] = df_exon['end'].astype(int) - df_exon['start'].astype(int) - 1
    df_exon.loc[df_exon.strand == '-', 'rel_pos'] = df_exon['length'] - df_exon['rel_pos'] + 1
    # find longest exon and choose increase step
    # longest exon = 50294 nts, step = 10**-5
    dec = 5
    df_exon['rel_pos_start'] = df_exon.rel_pos - 1
    df_exon['rel_pos_end'] = df_exon.rel_pos_start + 1
    df_exon['rel_pos_start'] /= df_exon['length']
    df_exon['rel_pos_end'] /= df_exon['length']

    df_exon = df_exon[['rel_pos_start', 'rel_pos_end', 'counts']].groupby(
        ['rel_pos_start', 'rel_pos_end']).counts.sum().reset_index()

    # normalize positions between 0 and 1
    step = 10 ** -dec
    pos = np.arange(0, 1 + step, step)
    zero = [0] * len(pos)
    init_count = [[pos[i], zero[i]] for i in range(len(pos))]
    df_init = pd.DataFrame(init_count, columns=['rel_pos', 'counts'])
    for idx, row in df_exon.iterrows():
        df_init.loc[(df_init.rel_pos >= row['rel_pos_start']) &
                    (df_init.rel_pos < row['rel_pos_end']), 'counts'] += row['counts']
    df_init['counts'] = df_init['counts'] / nb_exons * norm_factor
    return df_init


def draw_subplots(df, ax):
    ax.step(x=df.rel_pos, y=df.counts, where='post')
    ax.yaxis.set_tick_params(labelbottom=True)
    return ax


def plot_coverage(df_up, df_exon, df_down, title, outpath, pct_int):
    fig, axes = plt.subplots(1, 3, sharey='all')
    fig.set_size_inches(18, 7)
    fig.suptitle(title)
    axes[0].set_xlabel('Nucleotides upstream of exon')
    axes[0].set_ylabel('Normalized counts\n(Proportion of interactions/nb features*10E6)')
    axes[1].set_xlabel('Relative position in exon')
    axes[2].set_xlabel('Nucleotides downstream of exon')
    axes[0] = draw_subplots(df_up, axes[0])
    axes[1] = draw_subplots(df_exon[df_exon.rel_pos != 1], axes[1])
    axes[2] = draw_subplots(df_down, axes[2])
    axes[0].set_ylim(bottom=0)
    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.2)
    try:
        snoname = title.split()[1]
    except IndexError:
        snoname = 'all'
    txt = '%.2f%% of %s\npredicted interactions\nfall into these categories' % (pct_int, snoname)
    plt.tight_layout(rect=[0, 0.03, 0.98, 0.9])
    axes[2].text(1, 1.02, txt, transform=axes[2].transAxes, fontsize=10,
        verticalalignment='bottom', horizontalalignment='right')
    fig.savefig(os.path.join(outpath, title.split(' ')[0] + '.prop.svg'))
    plt.close()


def calc_sno_coverage(df, df_int, region):
    # bedtools coverage -d -counts
    # df feature = a, df interactions = b
    if region == 'exon':
        names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'rel_pos', 'counts']
    else:
        names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'block_length', 'rel_pos', 'counts']
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    bed_ft = BedTool.from_dataframe(df)
    bed_int = BedTool.from_dataframe(df_int)
    # only compute coverage of feature having at least one interaction
    bed_ft = bed_ft.intersect(bed_int, wa=True, s=True, u=True)
    cov = bed_ft.coverage(b=bed_int, s=True, d=True)
    df_cov = cov.to_dataframe(names=names, dtype={'chrom': str})
    if region != 'exon':
        df_cov = df_cov.drop(['block_length'], axis=1)
    return df_cov


def merge_df(df_all, df_to_add):
    df_all = pd.concat([df_all, df_to_add], ignore_index=True)
    df_all = df_all.groupby('rel_pos').counts.sum().reset_index()
    return df_all


def compute_sno(args):
    snoid, df_up, df_exon, df_down, inpath, len_dict, norm_factor, nb_windows, nb_exons, outpath, snoname = args
    int_file = os.path.join(inpath, 'pred_%s.98_3.gene.htrri.tsv' % (snoid))
    df_int = pd.read_csv(int_file, sep='\t',
                         names=['seqname', 'target_window_start', 'target_window_end', 'snoid', 'count',
                                'strand', 'sno_window_start', 'sno_window_end', 'mean_score', 'min_score',
                                'max_score', 'target'],
                         dtype={'seqname': str})
    df_int = df_int[df_int.seqname != 'seqname']  # remove header if present
    df_int[['target_window_start', 'target_window_end']] = df_int[['target_window_start',
                                                                   'target_window_end']].astype(int)
    df_int = df_int.drop(['target'], axis=1)
    df_int = df_int.drop_duplicates()
    df_int = df_int.drop(['sno_window_start', 'sno_window_end', 'mean_score', 'min_score',
                          'max_score'], axis=1)
    df_int[['count']] = df_int[['count']].astype(float)
    df_int = df_int[df_int['count'] >= nb_windows]
    df_int['target_window_end'] += 1
    tot_int = len(df_int)
    if  os.path.exists(os.path.join(outpath, snoid + '_down.csv')):
        sno_cov_up = pd.read_csv(os.path.join(outpath, snoid + '_up.csv'))
        sno_cov_exon = pd.read_csv(os.path.join(outpath, snoid + '_exon.csv'))
        sno_cov_down = pd.read_csv(os.path.join(outpath, snoid + '_down.csv'))
    else:
        cov_up = calc_sno_coverage(df_up, df_int, 'up')
        sno_cov_up = prep_up(cov_up, len_dict, norm_factor)
        cov_exon = calc_sno_coverage(df_exon, df_int, 'exon')
        # exon: normalize between 0 and 1, normalize counts like up and down
        sno_cov_exon = prep_exon(cov_exon, nb_exons, norm_factor)
        cov_down = calc_sno_coverage(df_down, df_int, 'down')
        sno_cov_down = prep_down(cov_down, len_dict, norm_factor)

        sno_cov_up.to_csv(os.path.join(outpath, snoid + '_up.csv'), index=False)
        sno_cov_exon.to_csv(os.path.join(outpath, snoid + '_exon.csv'), index=False)
        sno_cov_down.to_csv(os.path.join(outpath, snoid + '_down.csv'), index=False)

    bed_int = BedTool.from_dataframe(df_int[['seqname', 'target_window_start', 'target_window_end',
                                             'snoid', 'count', 'strand']])
    df_all = pd.concat([df_up, df_exon, df_down], sort=True)[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    df_all[['start', 'end']] = df_all[['start', 'end']].astype(int)
    bed_all = BedTool.from_dataframe(df_all)
    bed_res = bed_int.intersect(bed_all, wa=True, s=True, u=True)
    nb_int = bed_res.count()
    pct_int = nb_int / tot_int * 100

    plot_coverage(sno_cov_up, sno_cov_exon, sno_cov_down, snoid + ' ' + snoname, outpath, pct_int)
    return sno_cov_up , sno_cov_exon, sno_cov_down, tot_int, nb_int


def calc_coverage(df_up, df_exon, df_down, inpath, outdir, nb_windows, snolist, nb_threads, sno_dict):
    norm_factor = 10E6
    outpath = os.path.join(outdir, 'exon_flanking_region')
    os.makedirs(outpath, exist_ok=True)

    len_dict = {}
    len_dict['up'] = df_up.groupby('length').name.count().to_dict()
    len_dict['down'] = df_down.groupby('length').name.count().to_dict()

    all_cov_up = pd.DataFrame(columns=['rel_pos', 'counts'])
    all_cov_down = pd.DataFrame(columns=['rel_pos', 'counts'])
    all_cov_exon =  pd.DataFrame(columns=['rel_pos', 'counts'])
    all_tot = 0
    all_int = 0

    nb_exons = len(df_exon)
    with Pool(nb_threads) as p:
        res = p.map(compute_sno, [[snoid, df_up, df_exon, df_down, inpath, len_dict, norm_factor,
                                   nb_windows, nb_exons, outpath, sno_dict[snoid]] for snoid in snolist])
    for res_i in res:
        all_cov_up = merge_df(all_cov_up, res_i[0])
        all_cov_exon = merge_df(all_cov_exon, res_i[1])
        all_cov_down = merge_df(all_cov_down, res_i[2])
        all_tot += res_i[3]
        all_int += res_i[4]
    pct_int = all_int / all_tot * 100
    plot_coverage(all_cov_up, all_cov_exon, all_cov_down, 'All', outpath, pct_int)


def bed_to_df(bed, side):
    """
    create bed of up- or down- stream exon regions.
    :param bed: bed of exon positions
    :param side: up or down
    :return: dataframe with region up/down stream of exons
    """
    names = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    names += [side + '_'  + c for c in names]
    names += ['length']
    to_keep = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'length']
    df = bed.to_dataframe(names=names, dtype={'chrom': str})
    df = df[to_keep]

    # keep min(100 nts, 1/2 intron) to avoid having same part of an intron ni 5' and 3' of different exons
    df['length'] = df['length'].abs()
    df = df[df.length != 0] # remove 1st or last exon depending if up/down
    df.loc[df.length < 200, 'new_length'] = df.length / 2
    df.loc[df.length >= 200, 'new_length'] = 100
    df['length'] = df['new_length']
    df = df.drop(['new_length'], axis=1)
    df['length'] = df['length'].astype(int)
    if side == 'up':
        df.loc[df.strand == '+', 'end'] = df.start
        df.loc[df.strand == '+', 'start'] -= df.length
        df.loc[df.strand == '-', 'start'] = df.end
        df.loc[df.strand == '-', 'end'] += df.length
    elif side == 'down':
        df.loc[df.strand == '+', 'start'] = df.end
        df.loc[df.strand == '+', 'end'] = df.start + df.length
        df.loc[df.strand == '-', 'end'] = df.start
        df.loc[df.strand == '-', 'start'] = df.end - df.length
    return df


def main():
    gtf_file = sys.argv[1] # full human gtf file in csv format
    inpath = os.path.abspath(sys.argv[2]) # path to predicted interaction files
    outdir = os.path.abspath(sys.argv[3]) # path to output directory
    nb_windows = int(sys.argv[4]) # minimal number of consecutive windows
    snofile = sys.argv[5] #csv file with chromosome and sno_id
    nb_threads = int(sys.argv[6]) # nb of threads to use in parallel

    with open(snofile, 'r') as f:
        snolist = f.readlines()
    snolist = [i.strip().split(',')[1] for i in snolist]
    for sno in ['NR_003072', 'NR_003073', 'ENSG00000252284']: # SNORD91A, SNORD91B, SNORA108
        try:
            snolist.remove(sno)
        except ValueError:
            pass

    # read_gtf
    df_gtf = pd.read_csv(gtf_file, sep='\t',dtype={'seqname': str})
    chromo_dict = df_gtf[df_gtf.feature == 'gene'].set_index('gene_id').seqname.to_dict()
    df_sno = df_gtf[(df_gtf.feature == 'gene') & (df_gtf.gene_biotype == 'snoRNA')][['gene_id', 'gene_name']]
    sno_dict = df_sno.set_index('gene_id').gene_name.to_dict()

    # make exon bed
    df_exon = df_gtf[(df_gtf.feature == 'exon') & (df_gtf.gene_biotype == 'protein_coding')].copy(deep=True)
    df_exon['score'] = 0
    df_exon['name'] = df_exon.transcript_id + '|' + df_exon.exon_number.astype(str)
    df_exon = df_exon.drop_duplicates(subset=['gene_id', 'start', 'end', 'strand'])
    df_exon = df_exon[['gene_id', 'start', 'end', 'name', 'score', 'strand']]
    df_exon = df_exon.sort_values(['gene_id', 'start', 'end'])
    bed_exon = BedTool.from_dataframe(df_exon)

    # use bedtools closest to get closest up/downstream exons for each exon
    # bedtools closest correctly handles strand, upstream = 5' , downstream = 3'
    # upstream length is negative, downstream is positive
    bed_up = bed_exon.closest(bed_exon, s=True, D='a', t='first', fu=True)
    bed_down = bed_exon.closest(bed_exon, s=True, D='a', t='first', fd=True)

    df_up = bed_to_df(bed_up, 'up')
    df_down = bed_to_df(bed_down, 'down')

    # remove monoexonic genes
    id_list = df_up.name.unique().tolist()
    id_list.extend(df_down.name.unique().tolist())
    df_exon = df_exon[df_exon.name.isin(id_list)]

    # map chromo to gene_id
    df_up['chrom'] = df_up.chrom.map(chromo_dict)
    df_down['chrom'] = df_down.chrom.map(chromo_dict)
    df_exon['chrom'] = df_exon.gene_id.map(chromo_dict)
    df_exon = df_exon[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    calc_coverage(df_up, df_exon, df_down, inpath, outdir, nb_windows, snolist, nb_threads, sno_dict)


if __name__ == '__main__':
    main()