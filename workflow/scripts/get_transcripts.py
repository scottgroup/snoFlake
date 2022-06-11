#!/usr/bin/env python

"""
Get transcripts from transcriptome annotation file
[1] all: Get all transcripts
[2] rbp: Get pre-mRNA transcripts for RBPs
[3] sno: Get snoRNA transcripts
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
    return pd.read_csv(node_list, sep='\t')

def rbp_transcript(df_gtf, df_rbp): # get transcripts for RBPs only
    df = pd.DataFrame(columns=df_gtf.columns)
    for i in range(len(df_rbp)):
        curr_rbp = df_rbp.loc[i,'name']
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
        curr_sno = df_sno.loc[i,'id']
        df_temp = df_gtf[df_gtf["gene_id"]==curr_sno]
        df = pd.concat([df,df_temp], ignore_index=True)
    df = df[['seqname','start','end','gene_id','score','strand']]
    df.drop_duplicates(inplace=True,ignore_index=True)
    return df

def main():
    option = snakemake.params[0] # choose from: all, rbp, sno
    annotation = snakemake.input[0] # genome annotation in gtf format
    rbp = snakemake.input[1] # RBP list in tsv 
    sno = snakemake.input[2] # snoRNA list in tsv
    out = snakemake.output[0] # output file

    transcripts = read_annotation(annotation)

    if option == "all": # all transcripts
        all_transcripts(transcripts).to_csv(out,sep='\t',index=None,header=False)
    elif option == "rbp": # RBP transcripts
        rbp_lst = read_node_list(rbp)
        rbp_transcript(transcripts, rbp_lst).to_csv(out,sep='\t',index=None,header=False)
    elif option == "sno": # snoRNA transcripts
        sno_lst = read_node_list(sno)
        sno_transcript(transcripts, sno_lst).to_csv(out,sep='\t',index=None,header=False)
    else: # not a valid argument
        print("Transcript argument type not valid. Choose between ALL, RBP and SNO in lower case.")

if __name__ == '__main__':
    main()