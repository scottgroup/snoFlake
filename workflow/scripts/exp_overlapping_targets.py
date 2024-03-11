#!/usr/bin/env python


import pandas as pd
from pybedtools import BedTool
import os


def bedtools_intersect(interactions_df,ENCODE_dir,sno_bed,exp_pc_bed,outdir):
    """
    Bedtools intersect ENCODE RBP interactions with snoRNA binding interactions to get their overlapping targets.
    """
    for i in range(len(interactions_df)):
        RBP = interactions_df.loc[i,'RBP']
        RBP_bed = BedTool(os.path.join(ENCODE_dir, RBP+'.bed'))
        intersect_bed = sno_bed.intersect(RBP_bed, s=True)
        exp_overlapping_targets_bed = intersect_bed.intersect(exp_pc_bed,s=True,wo=True)
        print(exp_overlapping_targets_bed)
        column_names = ['intersect_Chromosome','intersect_Start','intersect_End','snoRNA','intersect_Score','intersect_Strand',
                                'pc_Chromosome','pc_Start','pc_End','pc_gene_id','pc_Score','pc_Strand']
        exp_overlapping_targets_df = exp_overlapping_targets_bed.to_dataframe(header=None,names=column_names)
        #exp_overlapping_targets_df = exp_overlapping_targets_df[['intersect_Chromosome','intersect_Start','intersect_End','pc_gene_id','intersect_Score','intersect_Strand']]
        exp_overlapping_targets_df.to_csv(os.path.join(outdir, RBP+'.tsv'),sep='\t',index=None)
        break
    return 

def main():
    interactions = snakemake.input.sno_RBP_overlap # all snoRNA-RBP target overlap interactions for one snoRNA
    p_val_thres = snakemake.params.p_val_thres # p-value threshold for significant snoRNA-RBP target overlap interactions
    outdir = snakemake.params.outdir
    exp_pc_genes = snakemake.params.exp_pc_genes # list of expressed protein-coding genes
    sno_interaction = snakemake.input.sno_interaction
    ENCODE_dir = snakemake.params.preprocessed_ENCODE # preprocessed ENCODE RBP interactions

    # Filter significant interactions
    interactions_df = pd.read_csv(interactions,delim_whitespace=True,names=['snoRNA','RBP','num1','num2'])
    interactions_df = interactions_df[(interactions_df['num1']==0) & (interactions_df['num2']==p_val_thres)]

    # Bedtools intersect each snoRNA-RBP pair
    sno_bed = BedTool(sno_interaction)
    exp_pc_df = pd.read_csv(exp_pc_genes,sep='\t')
    exp_pc_df['score'] = 100
    exp_pc_df = exp_pc_df[['Chromosome','Start','End','gene_id','score','Strand']]
    exp_pc_bed = BedTool.from_dataframe(exp_pc_df)

    bedtools_intersect(interactions_df,ENCODE_dir,sno_bed,exp_pc_bed,outdir)

    sum_file = os.path.join(outdir,'summary.txt')
    f = open(sum_file, "w")
    f.write("Overlapping targets obtained!")
    f.close()


if __name__ == '__main__':
    main()