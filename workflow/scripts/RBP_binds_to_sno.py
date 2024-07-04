#!/usr/bin/env python


import pandas as pd
import os
from pybedtools import BedTool
import sys


def bedtools_intersect(ENCODE_dir, sno_bed, fasta):
    """
    Bedtools intersect ENCODE RBP interactions with snoRNA location.
    """
    RBP_list = os.listdir(ENCODE_dir) # ENCODE RBP interaction directory

    final_interactions_df = pd.DataFrame(columns=['source','target','raw_score'])

    for RBP_file in RBP_list:
        RBP_bed = BedTool(os.path.join(ENCODE_dir, RBP_file))
        intersect_bed = RBP_bed.intersect(sno_bed, s=True, wb=True)
        column_names = ['RBP_Chromosome','RBP_Start','RBP_End','RBP_gene_name','RBP_Score','RBP_Strand',
                                'sno_Chromosome','sno_Start','sno_End','sno_gene_id','sno_Score','sno_Strand']
        intersect_df = intersect_bed.to_dataframe(header=None,names=column_names)

        if len(intersect_df)>0:
            # Get snoRNA sequence of the region bound by RBP
            intersect_df['RBP-sno'] = intersect_df['RBP_gene_name'] + intersect_df['sno_gene_id']
            intersect_pos_bed = BedTool.from_dataframe(intersect_df[['RBP_Chromosome','RBP_Start','RBP_End','RBP-sno','RBP_Score','RBP_Strand']])
            intersect_seq = intersect_pos_bed.sequence(fi=fasta,s=True,tab=True,name=True)
            print(open(intersect_seq.seqfn).read())

            intersect_df = intersect_df.rename(columns={"RBP_gene_name": "source","sno_gene_id" : "target", "RBP_Score" : "ENCODE_log_pval"})
            intersect_df = intersect_df[['source','target','ENCODE_log_pval']]
            intersect_df = intersect_df.groupby(by=['source','target'], as_index = False).min()
            final_interactions_df = pd.concat([final_interactions_df,intersect_df],ignore_index=True)
                
    final_interactions_df['interaction'] = 'RBP_binds_to_snoRNA'
    final_interactions_df.drop_duplicates(inplace=True, ignore_index=True)
    return final_interactions_df


def main():
    ENCODE_dir = sys.argv[1] # path to preprocessed ENCODE datasets
    sno_annot = sys.argv[2] # snoRNA annotation
    fasta = sys.argv[3]
    outfile = sys.argv[4]

    sno_annot_df = pd.read_csv(sno_annot,sep='\t',usecols=['Chromosome','Start','End','gene_id','Strand'])
    sno_annot_df['Score'] = 3
    sno_annot_df = sno_annot_df[['Chromosome','Start','End','gene_id','Score','Strand']]
    sno_bed = BedTool.from_dataframe(sno_annot_df) # snoRNA location in bed format

    df_interaction = bedtools_intersect(ENCODE_dir,sno_bed,fasta)

    df_interaction = df_interaction[['source','target','ENCODE_log_pval','interaction']]

    df_interaction.to_csv(outfile, sep='\t',index=False)


if __name__ == '__main__':
    main()
