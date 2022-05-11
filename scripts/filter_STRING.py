#!/usr/bin/env python

""" Filter STRING interactions by only keeping RBP-RBP pairs specified in the RBP list """

import pandas as pd
import sys

def rbp_list(file): # 2 columns (protein_id, name)
    return pd.read_csv(file,sep='\t')

def string_df(file): # 7 columns (item_id_a, item_id_b, mode, action, is_directional, a_is_acting, score)
    return pd.read_csv(file,sep='\t')

def select_interactions(rbp_file, string_file):
    # create df
    df_rbp = rbp_list(rbp_file)
    df_string = string_df(string_file)
    df_interactions = pd.DataFrame(columns=df_string.columns) # store final interactions

    # extract interactions for each RBP
    for i in range(len(df_rbp)):
        rbp = df_rbp.loc[i,'protein_id'] # current RBP
        temp = df_string[df_string.loc[:, 'item_id_a']==rbp].reset_index(drop=True)
        for j in range(len(temp)):
            if temp.loc[j,'item_id_b'] in df_rbp.loc[:, 'protein_id'].values:
                name = df_rbp[df_rbp.loc[:,'protein_id']==temp.loc[j,'item_id_b']].reset_index(drop=True).loc[0,'name'] # want RBP name not id
                #################### NEED TO FIX FROM HERE ##############
                df_new = pd.DataFrame(data = {'RBP1':[df_rbp.loc[i,'name']],'RBP2':[name]}) 
                df_interactions = pd.concat([df_interactions, df_new],ignore_index=True)
    # drop duplicate rows
    df_interactions.drop_duplicates(ignore_index=True,inplace=True)
    # check for reverse cases?
    return df_interactions

def main():
    rbp = sys.argv[1] # rbp list with protein id
    string = sys.argv[2] # unzipped STRING interaction file
    out = sys.argv[3] # output file

    select_interactions(rbp,string).to_csv(out, sep='\t',index=False)

if __name__ == '__main__':
    main()