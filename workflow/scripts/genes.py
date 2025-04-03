#!/usr/bin/env python

import pandas as pd
import pyranges as pr


"""
Extract all protein-coding genes from gtf and save as a BED file.
"""


def main():

    df_gtf = pr.read_gtf(snakemake.input.gtf).df
    df_gene = df_gtf[(df_gtf['Feature'] == 'gene') & (df_gtf['gene_biotype'] == 'protein_coding')]
    df_gene['Score'] = 3
    df_gene = df_gene[['Chromosome','Start','End','gene_id','Score','Strand']]
    df_gene.to_csv(snakemake.output[0],sep='\t',index=None,header=None)


if __name__ == '__main__':
    main()