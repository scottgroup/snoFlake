#!/usr/bin/env python


import pandas as pd
#from pybedtools import BedTool, helpers
#import pyranges as pr
import sys
import os


def compare_df(df1,df2):
    df_merged = df1.merge(df2, how='left', indicator=True)
    df1_uniq = df_merged[df_merged['_merge']=='left_only']
    df1_uniq.drop(columns=['_merge'],inplace=True)
    return df1_uniq


def main():
    files_list = snakemake.input # in order: [lenient, intermediate, stringent]
    outfile = snakemake.output[0]
    thres = snakemake.params[0]
    
    df_lenient = pd.read_csv(files_list[0],sep='\t')
    df_intermediate = pd.read_csv(files_list[1],sep='\t')
    df_stringent = pd.read_csv(files_list[2],sep='\t')

    # final df with weights
    df_final = pd.DataFrame(columns=['source','target','interaction','weight','p_val_thres'])
    
    df_stringent['weight'] = "stringent"
    df_stringent['p_val_thres'] = 1/thres[2]
    df_final = pd.concat([df_final,df_stringent],ignore_index=True)

    df_intermediate_uniq = compare_df(df_intermediate,df_stringent)
    df_intermediate_uniq['weight'] = "intermediate"
    df_intermediate_uniq['p_val_thres'] = 1/thres[1]
    df_final = pd.concat([df_final,df_intermediate_uniq],ignore_index=True)

    df_intermediate_stringent = pd.concat([df_stringent,df_intermediate],ignore_index=True).drop_duplicates(ignore_index=True)
    df_lenient_uniq = compare_df(df_lenient,df_intermediate_stringent)
    df_lenient_uniq['weight'] = "lenient"
    df_lenient_uniq['p_val_thres'] = 1/thres[0]
    df_final = pd.concat([df_final,df_lenient_uniq],ignore_index=True)

    df_final.sort_values(by=['source']).drop_duplicates(subset=['source','target','interaction']).to_csv(outfile,sep='\t',index=None)


if __name__ == '__main__':
    main()