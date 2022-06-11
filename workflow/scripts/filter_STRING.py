#!/usr/bin/env python

""" Filter STRING interactions by only keeping RBP-RBP pairs specified in the RBP list """

import pandas as pd
import sys

def rbp_list(file): # 2 columns (name, protein_id)
    return pd.read_csv(file,sep='\t',usecols=['name','protein_id'])

def string_df(file): # 7 columns (item_id_a, item_id_b, mode, action, is_directional, a_is_acting, score)
    return pd.read_csv(file,sep='\t')

def select_interactions(rbp_file, string_file):
    # create df
    df_rbp = rbp_list(rbp_file)
    df_string = string_df(string_file)
    df_interactions = pd.DataFrame(columns=['RBP1','RBP2','mode','score']) # store final interactions

    # extract interactions for each RBP
    for i in range(len(df_rbp)):
        rbp = df_rbp.loc[i,'protein_id'] # current RBP protein id
        temp = df_string[df_string.loc[:, 'item_id_a']==rbp].reset_index(drop=True)
        for j in range(len(temp)):
            if temp.loc[j,'item_id_b'] in df_rbp.loc[:, 'protein_id'].values:
                # want RBP names not protein id
                name1 = df_rbp[df_rbp.loc[:,'protein_id']==temp.loc[j,'item_id_a']].reset_index(drop=True).loc[0,'name']
                name2 = df_rbp[df_rbp.loc[:,'protein_id']==temp.loc[j,'item_id_b']].reset_index(drop=True).loc[0,'name']
                df_new = pd.DataFrame(data = {'RBP1':[name1],'RBP2':[name2],'mode':[temp.loc[j,'mode']],'score':[temp.loc[j,'score']]}) 
                df_interactions = pd.concat([df_interactions, df_new],ignore_index=True)
    return df_interactions

def remove_duplicate_edges(df):
    df.drop_duplicates(ignore_index=True,inplace=True)
    # remove duplicate edges(ignore edge direction) by storing interaction pairs as sets (unordered)
    interaction_list = []
    df_in_list = df.values.tolist()
    df_result = pd.DataFrame(columns=df.columns) # resulting df without duplicate edges
    for row in df_in_list:
        curr_set = set(row)
        if curr_set not in interaction_list:
            interaction_list.append(curr_set)
            curr_df = pd.DataFrame([row],columns=df.columns)
            df_result = pd.concat([df_result,curr_df],ignore_index=True)
    df_result['interaction'] = "STRING"
    return df_result

def main():
    rbp = snakemake.input[1] # rbp list with protein id
    string = snakemake.input[0] # unzipped STRING interaction file
    out = snakemake.output[0] # output file

    df_result = select_interactions(rbp,string)
    remove_duplicate_edges(df_result).to_csv(out, sep='\t',index=False)

if __name__ == '__main__':
    main()