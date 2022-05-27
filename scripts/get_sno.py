#!/usr/bin/python

import pandas as pd
import sys
from rank_sno import get_tpm

""" Create a list of snoRNAs from the annotation and keepy highly expressed ones. """

def sno_in_annotation(file):
    # extract snoRNAs from annotation
    df = pd.read_csv(file,sep='\t',usecols=['feature','gene_biotype','gene_name','gene_id'])
    df = df[df['feature']=='transcript']
    df = df[df['gene_biotype']=='snoRNA']
    df.drop(['feature', 'gene_biotype'], axis=1, inplace=True)
    df.drop_duplicates(inplace=True,ignore_index=True)
    return df

def filter_tpm(df,tpm):
    # filter snoRNAs found in annotation by TPM
    # get TPM for cell lines not in tissue
    df_tpm = pd.read_csv(tpm,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    # IF snoRNA TPM >= 10 in at least one cell line, THEN the snoRNA is kept in the list
    sno = pd.DataFrame(columns=['gene_id','gene_name'])
    for i in range(len(df)):
        _,_,max = get_tpm(df.loc[i,'gene_id'],df_tpm)
        if max >= 10:
            curr_sno = pd.DataFrame({'gene_id':[df.loc[i,'gene_id']],'gene_name':[df.loc[i,'gene_name']]})
            sno = pd.concat([sno,curr_sno],ignore_index=True)
    return sno

def filter_box(df,snodb):
    # snoDB in df
    df_snodb = pd.read_csv(snodb,sep='\t',usecols=['ensembl_id','refseq_id','gene_name','box_type'])

    # filter out snoRNAs with weird gene_ids
    df = df[~df['gene_id'].str.contains('cluster',na=False)].reset_index(drop=True)
    df = df[~df['gene_id'].str.contains('snoDB',na=False)].reset_index(drop=True)

    # remove box H/ACA snoRNAs, U3, U8,SNORD116
    for n in ['SNORA','U3','U8','SNORD116','SCA']:
        df = df[~df['gene_name'].str.contains(n,na=False)].reset_index(drop=True)

    # check box type for snoRNAs with names not starting with 'SNORD'
    cd = df[df['gene_name'].str.contains('SNORD',na=False)].reset_index(drop=True)
    others = df[~df['gene_name'].str.contains('SNORD',na=False)].reset_index(drop=True)
    # check box type from snoDB for others
    for i in range(len(others)):
        id = others.loc[i,'gene_id']
        info = df_snodb[df_snodb['ensembl_id'].str.contains(id,na=False)].reset_index(drop=True)
        if len(info) == 0:
            info = df_snodb[df_snodb['refseq_id']==id]
        if len(info) != 0 and info.loc[0,'box_type'] == 'C/D':
            curr_sno = pd.DataFrame({'gene_id':[id],'gene_name':[others.loc[i,'gene_name']]})
            cd = pd.concat([cd,curr_sno],ignore_index=True)
    return cd

def main():
    annotation = sys.argv[1] # annotation tsv file 
    tpm = sys.argv[2] # coco_tpm.tsv
    snodb = sys.argv[3] # snoDB.tsv
    outfile = sys.argv[4] # output file path

    initial_sno = sno_in_annotation(annotation)
    tpm_10 = filter_tpm(initial_sno,tpm)
    filter_box(tpm_10,snodb).to_csv(outfile,sep='\t',index=None)

if __name__ == '__main__':
    main()