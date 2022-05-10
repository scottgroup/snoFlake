#!/usr/bin/env python

import pandas as pd
import sys

""" Find snoRNAs that have RBPs as host genes in snoDB """

def load_db(file):
    # Load snoDB
    df = pd.read_csv(file,sep='\t',usecols=['ensembl_id','host_gene_name','host_biotype'])
    # Format snoDB by only keeping host genes that are protein coding
    return df[df['host_biotype']=='protein_coding'].dropna().reset_index(drop=True)

def filter_sno_host(df,snofile,rbpfile):
    result = pd.DataFrame(columns=['ensembl_id','host_gene_name'])
    df.drop(columns=['host_biotype'],inplace=True) # all protein_coding
    # keep only snoRNAs of interest
    sno_list = pd.read_csv(snofile, sep='\t').id.values.tolist()
    rbp_list = pd.read_csv(rbpfile, sep='\t').name1.values.tolist()
    for sno in sno_list:
        df_sno = df[df["ensembl_id"].str.contains(sno,na=False)].reset_index(drop=True)
        if len(df_sno)>0 and df_sno.loc[0,'host_gene_name'] in rbp_list:
            row = {'ensembl_id':sno,'host_gene_name':df_sno.loc[0,'host_gene_name']}
            result = result.append(row, ignore_index=True)
    # add interaction type column
    result['interaction'] = 'rbp_as_host_gene'
    return result


def main():
    snodb = sys.argv[1]
    snolist = sys.argv[2]
    rbplist = sys.argv[3]
    outfile = sys.argv[4]

    result = filter_sno_host(load_db(snodb),snolist,rbplist)
    result.to_csv(outfile, sep='\t',header=True,index=False)

if __name__ == '__main__':
    main()