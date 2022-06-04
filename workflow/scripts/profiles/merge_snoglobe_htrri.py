#!/usr/bin/env python

""" Merge snoGloBe predictions and HTRRI datasets for generating snoRNA and target profiles """

import pandas as pd
import os
import sys

def read_htrri(file):
    df = pd.read_csv(file, sep='\t', usecols=['chr1','start1','end1','strand1','chr2','start2','end2','strand2',
                                    'single_id1','single_id2','name1','name2','biotype1','biotype2'])
    # remove interactions where snoRNA interacts with itself
    for i in range(len(df)):
        if df.loc[i,'single_id1'] == df.loc[i,'single_id2']:
            df.drop(labels=i,axis=0,inplace=True)
    df.reset_index(drop=True,inplace=True)
    # add additional columns with constant value to fit snoglobe format
    df['count']=3
    df['mean_score']=0.99
    df['min_score']=0.99
    df['max_score']=0.99
    return df

def read_snodb(file,id): # get snoRNA start & end position in the genome from snoDB
    df = pd.read_csv(file, sep='\t', usecols=['ensembl_id','chr','start','end','length'])
    df_sno = df[df["ensembl_id"]==id].reset_index(drop=True)
    if len(df_sno) == 0: # exact ensembl id match not found
        df_sno = df[df["ensembl_id"].str.contains(id,na=False)].reset_index(drop=True)
    sno_start = df_sno.loc[0,'start']
    sno_end = df_sno.loc[0,'end']
    sno_length = df_sno.loc[0,'length']
    return sno_start, sno_end, sno_length

def compute_sno_window(df,id,snodb):
    sno_start,sno_end,sno_length = read_snodb(snodb,id)
    result = pd.DataFrame(columns=['sno_window_start','sno_window_end'])
    window_start = 0
    window_end = 0
    for i in range(len(df)):
        if int(df.iloc[i,0]) <= sno_start: 
            window_start = 1
        else:
            window_start = df.iloc[i,0] - sno_start + 1
        if int(df.iloc[i,1]) >= sno_end:
            window_end = sno_length
        else:
            window_end = df.iloc[i,1] - df.iloc[i,0] + window_start
        curr = pd.DataFrame([[window_start,window_end]], columns=['sno_window_start','sno_window_end'])
        result = pd.concat([result,curr],ignore_index=True)
    return result

def format_htrri(df,sno_id, snodb):
    temp1 = pd.DataFrame(columns=['seqname', 'target_window_start', 'target_window_end',
                                    'snoid', 'count', 'strand', 'sno_window_start',
                                    'sno_window_end', 'mean_score', 'min_score', 'max_score', 'target'])
    temp2 = pd.DataFrame(columns=['seqname', 'target_window_start', 'target_window_end',
                                    'snoid', 'count', 'strand', 'sno_window_start',
                                    'sno_window_end', 'mean_score', 'min_score', 'max_score', 'target'])
    # get all HTRRI interactions for the current sno_id in single_id1 column
    sno1 = df[df["single_id1"]== sno_id].reset_index(drop=True)
    if len(sno1)>1:
        temp1['seqname'] = "chr"+sno1['chr2']
        temp1['target_window_start'] = sno1['start2']
        temp1['target_window_end'] = sno1['end2']
        temp1['snoid'] = sno1['single_id1']
        temp1['count'] = sno1['count']
        temp1['strand'] = sno1['strand2']
        # compute sno window
        df_computed = compute_sno_window(sno1[['start1','end1']],sno_id,snodb)
        temp1['sno_window_start'] = df_computed['sno_window_start']
        temp1['sno_window_end'] = df_computed['sno_window_end']
        temp1['mean_score'] = sno1['mean_score']
        temp1['min_score'] = sno1['min_score']
        temp1['max_score'] = sno1['max_score']
        temp1['target'] = sno1['single_id2']
    # get all HTRRI interactions for the current sno_id in single_id2 column
    sno2 = df[df["single_id2"]== sno_id].reset_index(drop=True)
    if len(sno2)>1:
        temp2['seqname'] = "chr"+sno2['chr1']
        temp2['target_window_start'] = sno2['start1']
        temp2['target_window_end'] = sno2['end1']
        temp2['snoid'] = sno2['single_id2']
        temp2['count'] = sno2['count']
        temp2['strand'] = sno2['strand1']
        # compute sno window
        df_computed = compute_sno_window(sno1[['start2','end2']],sno_id,snodb)
        temp2['sno_window_start'] = df_computed['sno_window_start']
        temp2['sno_window_end'] = df_computed['sno_window_end']
        temp2['mean_score'] = sno2['mean_score']
        temp2['min_score'] = sno2['min_score']
        temp2['max_score'] = sno2['max_score']
        temp2['target'] = sno2['single_id1']
    df_formatted = pd.concat([temp1,temp2],ignore_index=True)
    return df_formatted

def read_snoglobe(file):
    df = pd.read_csv(file, sep='\t', names=['seqname', 'target_window_start', 'target_window_end',
                                    'snoid', 'count', 'strand', 'sno_window_start',
                                    'sno_window_end', 'mean_score', 'min_score', 'max_score', 'target'])
    return df

def merge_interactions(snoglobe_file, df_htrri, snoid, snodb):
    df_snoglobe = read_snoglobe(snoglobe_file)
    df_htrri_formatted = format_htrri(df_htrri, snoid, snodb)
    return pd.concat([df_snoglobe, df_htrri_formatted], ignore_index=True)

def main():
    snoglobe = sys.argv[1] # path to snoglobe predictions directory
    htrri = sys.argv[2] # HTRRI file
    sno_file = sys.argv[3] # csv file containing chromosome, snoRNA ENSEMBL ID, snoRNA name
    snodb = sys.argv[4] # snoDB file
    outdir = sys.argv[5] # output directory to store merged snoRNA interaction files
    
    df_htrri = read_htrri(htrri)

    sno_list = pd.read_csv(sno_file,sep=',',names=['chromosome','ensembl_id','name']).ensembl_id.tolist()
    for sno in sno_list:
        print(sno)
        prediction_file = os.path.join(snoglobe,"pred_"+sno+".98_3.gene.tsv")
        out_file = os.path.join(outdir,"pred_"+sno+".98_3.gene.htrri.tsv")
        merge_interactions(prediction_file, df_htrri,sno,snodb).to_csv(out_file,sep='\t',index=None,header=None)

if __name__ == '__main__':
    main()