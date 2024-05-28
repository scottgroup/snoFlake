#!/usr/bin/env python


import pandas as pd
import numpy as np
from composite_score_ENCODE import min_max_normalize


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
    print('Significance threshold: ' + str(df.loc[len(df)-1,'right']))
    return df


def composite_score(df):
    """
    Compute a composite score that considers both the p-value and enrichment ratio.
    composite score = -log10(p-value) * enrichment ratio
    """
    # -log10 transform p-values
    df['log_right'] = -np.log10(df['right'])
    df['raw_composite_score'] = df['log_right'] * df['ratio']

    # separate df by p-value
    df_pval_zero = df[df['right']==0]
    df_pval_nonzero = df[df['right']>0]

    # normalize composite score
    df_pval_nonzero['normalized_composite_score'] = min_max_normalize(df_pval_nonzero['raw_composite_score'])

    # if p-value is 0, normalized_composite_score = 1
    p_val_zero = df_pval_zero.index.values.tolist()
    for i in p_val_zero:
        df_pval_zero.at[i, 'normalized_composite_score'] = 1

    df_merged = pd.concat([df_pval_zero,df_pval_nonzero],ignore_index=True)

    return df_merged


def format_df(df):
    """
    Format dataframe to match Cytoscape input style.
    """
    df['interaction'] = "sno_RBP_overlap"
    df = df[['snoRNA','RBP','normalized_composite_score','interaction']]
    df.rename(columns={'snoRNA': 'source', 'RBP': 'target', 'normalized_composite_score':'normalized_score'}, inplace=True)
    return df


def main():
    
    # sno_RBP_overlap df
    df = pd.read_csv(snakemake.params.merged,sep='\t',names=['left','right','two-tail','ratio','snoRNA','RBP'])
    df = df[['snoRNA','RBP','left','right','two-tail','ratio']]
    # desired false discovery rate
    fdr = snakemake.params.fdr

    df = format_df(composite_score(bh(df,fdr)))
    
    df.to_csv(snakemake.output.bh,sep='\t',index=None)


if __name__ == '__main__':
    main()