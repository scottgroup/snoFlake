#!/usr/bin/env python


import pandas as pd
import os
from pybedtools import BedTool


def bedtools_intersect(ENCODE_dir, sno_bed):
    """
    Bedtools intersect ENCODE RBP interactions with snoRNA location.
    """
    RBP_list = os.listdir(ENCODE_dir) # ENCODE RBP interaction directory

    final_interactions_df = pd.DataFrame(columns=['source','target'])

    for RBP_file in RBP_list:
        RBP_bed = BedTool(os.path.join(ENCODE_dir, RBP_file))
        intersect_bed = RBP_bed.intersect(sno_bed, s=True, wb=True)
        column_names = ['RBP_Chromosome','RBP_Start','RBP_End','RBP_gene_name','RBP_Score','RBP_Strand',
                                'sno_Chromosome','sno_Start','sno_End','sno_gene_id','sno_Score','sno_Strand']
        intersect_df = intersect_bed.to_dataframe(header=None,names=column_names)

        if len(intersect_df)>0:
            intersect_df = intersect_df.rename(columns={"RBP_gene_name": "source","sno_gene_id" : "target"})
            intersect_df = intersect_df[['source','target']]
            final_interactions_df = pd.concat([final_interactions_df,intersect_df],ignore_index=True)

            # get sequence on the snoRNA CHANGE MEEEEEEEEEEEEEEEEEEeeeeeeee ########
            intersect_bed = RBP_bed.intersect(sno_bed, s=True)
            intersect_df = intersect_bed.to_dataframe()
            print(intersect_dfst)
                

    final_interactions_df['interaction'] = 'RBP_binds_to_snoRNA'
    final_interactions_df.drop_duplicates(inplace=True, ignore_index=True)
    return final_interactions_df

def main():
    ENCODE_dir = snakemake.params.preprocessed_ENCODE # path to preprocessed ENCODE datasets
    sno_annot = snakemake.params.sno_annotation # snoRNA annotation
    outfile = snakemake.output[0]

    sno_annot_df = pd.read_csv(sno_annot,sep='\t',usecols=['Chromosome','Start','End','gene_id','Strand'])
    sno_annot_df['Score'] = 3
    sno_annot_df = sno_annot_df[['Chromosome','Start','End','gene_id','Score','Strand']]
    sno_bed = BedTool.from_dataframe(sno_annot_df) # snoRNA location in bed format

    bedtools_intersect(ENCODE_dir,sno_bed).to_csv(outfile, sep='\t',index=False)


if __name__ == '__main__':
    main()