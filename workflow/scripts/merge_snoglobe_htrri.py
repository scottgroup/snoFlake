#!/usr/bin/env python

import pandas as pd
from pybedtools import BedTool
import sys

"Merge snoGloBe predictions and HTRRI datasets."


def filter_htrri(file,sno,support_thres):

    df = pd.read_csv(file, sep='\t', usecols=['chr1','start1','end1','strand1','chr2','start2','end2','strand2',
                                    'gene_id1','gene_id2','gene_name1','gene_name2','gene_biotype1','gene_biotype2','support','exp'])
    # only want significant snoRNA-mRNA pairs
    support_thres = int(support_thres)
    df1 = df[(df["gene_id1"]==sno) & (df["gene_biotype2"]=="protein_coding") & (df["support"]>=support_thres)]
    df1 = df1[['chr2','start2','end2','gene_id1','support','strand2']]
    df1.rename(columns={'chr2':'chrom','start2':'start','end2':'end','gene_id1':'name','support':'score','strand2':'strand'},inplace=True)
    df2 = df[(df["gene_biotype1"]=="protein_coding") & (df["gene_id2"]==sno) & (df["support"]>=support_thres)]
    df2 = df2[['chr1','start1','end1','gene_id2','support','strand1']]
    df2.rename(columns={'chr1':'chrom','start1':'start','end1':'end','gene_id2':'name','support':'score','strand1':'strand'},inplace=True)
    df = pd.concat([df1,df2],ignore_index=True)
    df.reset_index(drop=True,inplace=True)
    return df


def filter_snoglobe(file,snoglobe_thres):
    snoglobe_thres = float(snoglobe_thres)
    df = pd.read_csv(file, sep='\t')
    df = df[df['min_score']>=snoglobe_thres]
    df = df[['target_chromo','target_window_start','target_window_end','sno_id','min_score','target_strand']]
    df.rename(columns={"target_chromo": "chrom", "target_window_start": "start", "target_window_end": "end", "sno_id": "name",
                        "min_score": "score", "target_strand": "strand"},inplace=True)
    return df


def main():
    snoglobe = sys.argv[1] # snoGloBe prediction file
    htrri = sys.argv[2] # HTRRI file
    sno = sys.argv[3] # current snoRNA ID
    outfile = sys.argv[4] # output file to store merged snoRNA interactions
    snoglobe_thres = sys.argv[5]
    support_thres = sys.argv[6]
    
    df_htrri = filter_htrri(htrri,sno,support_thres)
    df_snoglobe = filter_snoglobe(snoglobe,snoglobe_thres)

    df_merged = pd.concat([df_snoglobe,df_htrri],ignore_index=True)
    df_merged = BedTool.from_dataframe(df_merged).sort().merge(s=True,c=[4,5,6],o=['distinct','min','distinct']).to_dataframe()
    
    df_merged.to_csv(outfile,sep='\t',index=None,header=None)


if __name__ == '__main__':
    main()
