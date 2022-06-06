#!/usr/bin/env python

import pandas as pd

"Merge snoGloBe predictions and HTRRI datasets."


def read_htrri(file):
    df = pd.read_csv(file, sep='\t', usecols=['chr1','start1','end1','strand1','chr2','start2','end2','strand2',
                                    'single_id1','single_id2','name1','name2','biotype1','biotype2'])
    # only want snoRNA-RBP mRNA pairs
    df = df[(df["biotype1"]=="snoRNA") & (df["biotype2"]=="protein_coding")]
    df.reset_index(drop=True,inplace=True)

    # add additional columns with constant value to fit snoglobe format
    df['count']=3
    df['mean_score']=0.99
    df['min_score']=0.99
    df['max_score']=0.99
    return df

def read_snodb(file,id): # get snoRNA start & end position in the genome from snoDB
    df = pd.read_csv(file, sep='\t', usecols=['ensembl_id','refseq_id','chr','start','end','length'])
    df_sno = df[df["ensembl_id"].str.contains(id,na=False)].reset_index(drop=True)
    # id that starts with NR_
    if "NR_" in id:
        df_sno = df[df["refseq_id"].str.contains(id,na=False)].reset_index(drop=True)
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
    df_result = pd.DataFrame(columns=['seqname', 'target_window_start', 'target_window_end',
                                    'snoid', 'count', 'strand', 'sno_window_start',
                                    'sno_window_end', 'mean_score', 'min_score', 'max_score', 'target'])
    # get all HTRRI interactions for the current sno_id in single_id1 column
    sno = df[df["single_id1"]== sno_id].reset_index(drop=True)
    if len(sno)>1:
        df_result['seqname'] = "chr"+sno['chr2']
        df_result['target_window_start'] = sno['start2']
        df_result['target_window_end'] = sno['end2']
        df_result['snoid'] = sno['single_id1']
        df_result['count'] = sno['count']
        df_result['strand'] = sno['strand2']
        # compute sno window
        df_computed = compute_sno_window(sno[['start1','end1']],sno_id,snodb)
        df_result['sno_window_start'] = df_computed['sno_window_start']
        df_result['sno_window_end'] = df_computed['sno_window_end']
        df_result['mean_score'] = sno['mean_score']
        df_result['min_score'] = sno['min_score']
        df_result['max_score'] = sno['max_score']
        df_result['target'] = sno['single_id2']
    return df_result

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
    snoglobe = snakemake.input[0] # snoGloBe prediction file
    htrri = snakemake.params[0] # HTRRI file
    snodb = snakemake.params[1] # snoDB file
    sno = snakemake.params[2] # current snoRNA ID
    outfile = snakemake.output[0] # output file to store merged snoRNA interactions
    
    df_htrri = read_htrri(htrri)

    merge_interactions(snoglobe, df_htrri,sno,snodb).to_csv(outfile,sep='\t',index=None,header=None)

if __name__ == '__main__':
    main()