#!/usr/bin/env python


import pandas as pd
import numpy as np
import math


def bh(df,fdr):
    """
    Benjamini-Hochberg procedure to find the significance threshold.
    """
    # sort p-values from smallest to largest
    df['right'] = df['right'].astype(float)
    df = df.sort_values(by=['right']).reset_index(drop=True)

    # Benjamini-Hochberg critical value calculation
    m = len(df)
    bh_critical_vals = (np.arange(1, m + 1) / m) * float(fdr)
    df['bh_critical_vals'] = bh_critical_vals

    # find the largest i such that p(i) <= (i/m) * fdr
    df = df[df['right']<=df['bh_critical_vals']]

    # significance threshold
    sig_thres = df.loc[len(df)-1,'right'] # exact value
    print('Significance threshold: ' + str(sig_thres))
    sig_thres_exp = math.ceil(-np.log10(sig_thres)) # nearest largest exponent
    print('Significance threshold nearest exponent: -' + str(sig_thres_exp))

    df = df[df['right']<= 10 ** -sig_thres_exp]

    return df


def format_df(df):
    """
    Format dataframe to match Cytoscape input style.
    """
    df['interaction'] = "sno_RBP_overlap"
    df['sno_RBP_overlap_log_pval'] = -np.log10(df['right'])
    df = df[['snoRNA','RBP','sno_RBP_overlap_log_pval','interaction']]
    df.rename(columns={'snoRNA': 'source', 'RBP': 'target'}, inplace=True)
    return df


def main():
    
    # sno_RBP_overlap df
    df = pd.read_csv(snakemake.params.merged,sep='\t',names=['left','right','two-tail','ratio','snoRNA','RBP'])
    df = df[['snoRNA','RBP','left','right','two-tail','ratio']]
    # desired false discovery rate
    fdr = snakemake.params.fdr

    df = format_df(bh(df,fdr))
    df.to_csv(snakemake.output.bh,sep='\t',index=None)


if __name__ == '__main__':
    main()