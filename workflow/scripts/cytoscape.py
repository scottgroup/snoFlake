#!/usr/bin/env python

import pandas as pd
import sys


def node_weight(df):
    """
    Compute node size/weight based on TPM.
    """
    df['weight'] = 1
    df.loc[(df['max_TPM'] >= 10) & (df['max_TPM'] <100),'weight'] = 10
    df.loc[(df['max_TPM'] >= 100) & (df['max_TPM'] <1000),'weight'] = 100
    df.loc[df['max_TPM'] >= 1000,'weight'] = 1000
    return df


def edge_weight(df):
    """
    Compute edge thickness/weight based on p-values and scores.
    """
    df['weight'] = 1
    df.loc[(df['ENCODE_log_pval'] >= 10) & (df['ENCODE_log_pval'] <20),'weight'] = 10
    df.loc[(df['ENCODE_log_pval'] >= 20) & (df['ENCODE_log_pval'] <30),'weight'] = 100
    df.loc[df['ENCODE_log_pval'] >= 30,'weight'] = 1000
    df.loc[(df['sno_RBP_overlap_log_pval'] >= 10) & (df['sno_RBP_overlap_log_pval'] <20),'weight'] = 10
    df.loc[(df['sno_RBP_overlap_log_pval'] >= 20) & (df['sno_RBP_overlap_log_pval'] <30),'weight'] = 100
    df.loc[df['sno_RBP_overlap_log_pval'] >= 30,'weight'] = 1000
    df.loc[(df['STRING_score'] >= 925) & (df['STRING_score'] <950),'weight'] = 10
    df.loc[(df['STRING_score'] >= 950) & (df['STRING_score'] <975),'weight'] = 100
    df.loc[df['STRING_score'] >= 975,'weight'] = 1000
    return df


def main():

    # NODES
    sno = pd.read_csv(snakemake.input.snoRNA,sep='\t',usecols=['gene_id','gene_name','gene_biotype','min_TPM','max_TPM'])
    rbp = pd.read_csv(snakemake.input.RBP,sep='\t',usecols=['gene_name','gene_biotype','min_TPM','max_TPM'])
    rbp['gene_id'] = rbp['gene_name']
    rbp = rbp[['gene_id','gene_name','gene_biotype','min_TPM','max_TPM']]

    nodes = pd.concat([sno,rbp],ignore_index=True)
    nodes.gene_name.fillna(nodes.gene_id, inplace=True)
    nodes = node_weight(nodes)
    nodes.to_csv(snakemake.output.nodes,sep='\t',index=None)

    # EDGES
    rbp_binds_to_sno = pd.read_csv(snakemake.input.RBP_binds_to_sno[0],sep='\t')
    sno_rbp_ovlp = pd.read_csv(snakemake.input.sno_RBP_overlap,sep='\t')
    string = pd.read_csv(snakemake.input.STRING[0],sep='\t')
    edges = pd.concat([rbp_binds_to_sno,sno_rbp_ovlp,string],axis=0,ignore_index=True)
    edges = edges[['source','target','interaction','ENCODE_log_pval','sno_RBP_overlap_log_pval','STRING_score']]
    edges = edge_weight(edges)
    edges.to_csv(snakemake.output.edges,sep='\t',index=None)


if __name__ == '__main__':
    main()