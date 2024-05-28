#!/usr/bin/env python

import pandas as pd
import sys


def main():

    # NODES
    sno = pd.read_csv(snakemake.input.snoRNA,sep='\t',usecols=['gene_id','gene_name','gene_biotype','min_TPM','max_TPM'])
    rbp = pd.read_csv(snakemake.input.RBP,sep='\t',usecols=['gene_name','gene_biotype','min_TPM','max_TPM'])
    rbp['gene_id'] = rbp['gene_name']
    rbp = rbp[['gene_id','gene_name','gene_biotype','min_TPM','max_TPM']]

    nodes = pd.concat([sno,rbp],ignore_index=True)
    nodes.gene_name.fillna(nodes.gene_id, inplace=True)
    nodes.to_csv(snakemake.output.nodes,sep='\t',index=None)

    # EDGES
    rbp_binds_to_sno = pd.read_csv(snakemake.input.RBP_binds_to_sno[0],sep='\t')
    sno_rbp_ovlp = pd.read_csv(snakemake.input.sno_RBP_overlap,sep='\t')
    string = pd.read_csv(snakemake.input.STRING[0],sep='\t')
    edges = pd.concat([rbp_binds_to_sno,sno_rbp_ovlp,string],axis=0,ignore_index=True)
    edges.to_csv(snakemake.output.edges,sep='\t',index=None)


if __name__ == '__main__':
    main()