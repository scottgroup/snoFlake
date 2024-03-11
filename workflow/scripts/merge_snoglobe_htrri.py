#!/usr/bin/env python

import pandas as pd
import from pybedtools import BedTool

"Merge snoGloBe predictions and HTRRI datasets."


def filter_htrri(file,sno,support_thres):

    df = pd.read_csv(file, sep='\t', usecols=['chr1','start1','end1','strand1','chr2','start2','end2','strand2',
                                    'gene_id1','gene_id2','gene_name1','gene_name2','gene_biotype1','gene_biotype2','support','exp'])
    # only want significant snoRNA-mRNA pairs (support >= 5)
    df1 = df[(df["gene_id1"]==sno) & (df["gene_biotype2"]=="protein_coding") & (df["support"]>=support_thres)]
    df2 = df[(df["gene_biotype1"]=="protein_coding") & (df["gene_id2"]==sno) & (df["support"]>=support_thres)]
    df = pd.concat([df1,df2],ignore_index=True)
    df.reset_index(drop=True,inplace=True)

    return df


def filter_snoglobe(file,snoglobe_thres):
    df = pd.read_csv(file, sep='\t')
    df = df[df['min_score']>=snoglobe_thres]
    df = df[['target_chromo','target_window_start','target_window_end','sno_id','count','target_strand']]
    return df


def merge_interactions(df_snoglobe, df_htrri):
    snoglobe_bed = BedTool.from_dataframe(df_snoglobe)
    htrri_bed = BedTool.from_dataframe(df_htrri)
    merged_bed = snoglobe_bed.cat(htrri_bed,postmerge=False).sort().merge(s=True,c=[4,5,6],o=['distinct','max','distinct'])
    return merged_bed


def main():
    snoglobe = snakemake.input.snoglobe # snoGloBe prediction file
    htrri = snakemake.input.htrri # HTRRI file
    sno = snakemake.params.sno # current snoRNA ID
    outfile = snakemake.output[0] # output file to store merged snoRNA interactions
    snoglobe_thres = snakemake.params.snoglobe_thres
    support_thres = snakemake.params.support_thres
    
    df_htrri = filter_htrri(htrri,sno,support_thres)
    df_snoglobe = filter_snoglobe(snoglobe,snoglobe_thres)

    df_merged = pd.DataFrame(columns=['chrom', 'start', 'stop', 'name', 'score', 'strand'])

    if len(df_htrri)>0:
        df_merged = merge_interactions(df_snoglobe, df_htrri).to_dataframe()
    else:
        df_merged = df_snoglobe

    df_merged.to_csv(outfile,sep='\t',index=None,header=None)


if __name__ == '__main__':
    main()