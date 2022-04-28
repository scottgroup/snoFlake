#!/usr/bin/env python

""" Find snoRNA-snoRNA and snoRNA-RBP mRNA interactions from LIGR-seq, PARIS, SPLASH  datasets """

import pandas as pd
import sys

def read_htrri(file):
    df = pd.read_csv(file, sep='\t', usecols=['single_id1','single_id2','name1','name2','biotype1','biotype2'])
    # only want snoRNA-snoRNA and snoRNA-RBP mRNA pairs
    df = df[df["biotype1"]=="snoRNA"]
    df1 = df[df["biotype2"]=="snoRNA"]
    df2 = df[df["biotype2"]=="protein_coding"]
    return pd.concat([df1,df2],ignore_index=True)

def read_sno_rbp(sno_file, rbp_file):
    df_sno = pd.read_csv(sno_file, sep='\t')
    df_rbp = pd.read_csv(rbp_file, sep='\t')
    return df_sno, df_rbp

def get_interactions(df_htrri,df_sno,df_rbp):
    df_interactions = pd.DataFrame(columns=['RNA1','RNA2']) # store final interactions
    for i in range(len(df_sno)):
        sno = df_sno.iloc[i,0] # current sno
        temp = df_htrri[df_htrri['single_id1']==sno].reset_index(drop=True)
        for j in range(len(temp)):
            if temp.loc[j,'biotype2'] == "snoRNA": # RNA2 == snoRNA
                if temp.loc[j,'single_id2'] in df_sno.values:
                    df_new = pd.DataFrame(data = {'RNA1':[temp.loc[j,'single_id1']],'RNA2':[temp.loc[j,'single_id2']]})
                    df_interactions = pd.concat([df_interactions, df_new],ignore_index=True)
            else: # RNA2 == protein_coding
                if temp.loc[j,'name2'] in df_rbp.values:
                    df_new = pd.DataFrame(data = {'RNA1':[temp.loc[j,'single_id1']],'RNA2':[temp.loc[j,'name2']]})   
                    df_interactions = pd.concat([df_interactions, df_new],ignore_index=True)
    # add interaction type column
    df_interactions['interaction'] = 'htrri'
    # drop duplicate rows
    df_interactions.drop_duplicates(ignore_index=True,inplace=True)
    return df_interactions

def main():
    htrri = sys.argv[1]
    rbp = sys.argv[2]
    sno = sys.argv[3]
    out = sys.argv[4]

    df_htrri = read_htrri(htrri)
    df_sno, df_rbp = read_sno_rbp(sno, rbp)
    df_filtered = get_interactions(df_htrri,df_sno,df_rbp)
    df_filtered.to_csv(out, sep='\t',index=False)

if __name__ == '__main__':
    main()