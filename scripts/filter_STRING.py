#!/usr/bin/env python

""" Filter STRING interactions by only keeping RBP-RBP pairs specified in the RBP list """

import pandas as pd
import sys

def rbp_list(file): # 2 columns
    return pd.read_csv(file,sep='\t',header=None)

def string_df(file): # 2 columns
    return pd.read_csv(file,sep='\t',header=None)

def select_interactions(df_rbp, df_string):
    df_interactions = pd.DataFrame(columns=['RBP1','RBP2']) # store final interactions
    for i in range(len(df_rbp)):
        rbp = df_rbp.iloc[i,1] # current RBP
        temp = df_string[df_string.iloc[:, 0]==rbp].reset_index(drop=True)
        for j in range(len(temp)):
            if temp.iloc[j,1] in df_rbp.iloc[:, 1].values:
                name = df_rbp[df_rbp.iloc[:,1]==temp.iloc[j,1]].iloc[0,0] # want RBP name not id
                df_new = pd.DataFrame(data = {'RBP1':[df_rbp.iloc[i,0]],'RBP2':[name]})
                df_interactions = pd.concat([df_interactions, df_new],ignore_index=True)
    # add interaction type column
    df_interactions['interaction'] = 'STRING'
    # drop duplicate rows
    df_interactions.drop_duplicates(ignore_index=True,inplace=True)
    return df_interactions

def main():
    rbp = sys.argv[1]
    string = sys.argv[2]
    out = sys.argv[3]

    df_rbp = rbp_list(rbp)
    df_string = string_df(string)
    df_final = select_interactions(df_rbp,df_string)
    df_final.to_csv(out, sep='\t',index=False)

if __name__ == '__main__':
    main()