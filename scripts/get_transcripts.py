#!/usr/bin/env python

"""
Get transcripts from transcriptome annotation file
[1] Get all transcripts
[2] Get pre-mRNA transcripts for RBPs
[3] Get snoRNA transcripts
"""

import pandas as pd
pd.options.mode.chained_assignment = None
import sys

def read_annotation(gtf):
    df = pd.read_csv(gtf, sep='\t', usecols=['seqname','feature','start','end','gene_name','score',
        'strand','gene_id'])
    df = df[df['feature']=="transcript"] # keep transcripts only
    df = df.drop(columns=['feature'])
    return df

def read_node_list(node_list): # get RBP or snoRNA list
    return pd.read_csv(node_list, sep='\t',header=None)

def rbp_transcript(df_gtf, df_rbp): # get transcripts for RBPs only
    df = pd.DataFrame(columns=df_gtf.columns)
    for i in range(len(df_rbp)):
        curr_rbp = df_rbp.iloc[i,0]
        if curr_rbp == "TROVE2": # TROVE2 exists as RO60 in annotation (synonyms)
            curr_rbp = "RO60"
        df_temp = df_gtf[df_gtf["gene_name"]==curr_rbp]
        df = pd.concat([df,df_temp], ignore_index=True)

    # replace RO60 back to TROVE2
    df = df.replace("RO60","TROVE2")
    # change column order
    df = df[['seqname','start','end','gene_name','score','strand']]
    df.drop_duplicates(inplace=True,ignore_index=True)
    return df

def all_transcripts(df): # get all trancripts in annotation file
    df = df[['seqname','start','end','gene_name','score','strand']]
    df.drop_duplicates(inplace=True,ignore_index=True)
    return df

def sno_transcript(df_gtf, df_sno): # get transcripts for snoRNA only
    df = pd.DataFrame(columns=df_gtf.columns)
    for i in range(len(df_sno)):
        curr_sno = df_sno.iloc[i,0]
        df_temp = df_gtf[df_gtf["gene_id"]==curr_sno]
        df = pd.concat([df,df_temp], ignore_index=True)
    df = df[['seqname','start','end','gene_id','score','strand']]
    df.drop_duplicates(inplace=True,ignore_index=True)
    return df

def main():
    option = sys.argv[1]
    annotation = sys.argv[2]
    rbp = sys.argv[3]
    sno = sys.argv[4]
    out = sys.argv[5]

    transcripts = read_annotation(annotation)

    if option == str(1): # all transcripts
        df_result = all_transcripts(transcripts)
    elif option == str(2): # RBP transcripts
        rbp_lst = read_node_list(rbp)
        df_result = rbp_transcript(transcripts, rbp_lst)
    else: # option == 3 --> snoRNA transcripts
        sno_lst = read_node_list(sno)
        df_result = sno_transcript(transcripts, sno_lst)
    
    df_result.to_csv(out,sep='\t',index=None,header=False)

if __name__ == '__main__':
    main()