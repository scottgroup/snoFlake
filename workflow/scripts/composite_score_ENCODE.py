#!/usr/bin/env python

import pandas as pd
from pybedtools import BedTool
import sys


"""
Compute composite score for ENCODE eCLIP datasets using their p-values and signal values.
"""


def min_max_normalize(series):
    return (series - series.min()) / (series.max() - series.min())


def composite_score(df):
    """
    Compute composite score = -log10(p-value) * signal value
    """
    df['raw_composite_score'] = df['pValue'] * df['signalValue']
    return df


def main():
    cols = ['chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak']
    df = pd.read_csv(snakemake.input[0],sep='\t',names=cols)

    # composite score
    df = composite_score(df)
    df = df[['chrom','chromStart','chromEnd','name','raw_composite_score','strand']]

    # bedtools sort
    df = BedTool.from_dataframe(df).sort().to_dataframe()
    df.to_csv(snakemake.output[0],sep='\t',index=None,header=None)


if __name__ == '__main__':
    main()