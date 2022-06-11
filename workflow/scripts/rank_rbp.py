#!/usr/bin/python

from cmath import nan
from doctest import DocFileSuite
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np

""" Filter and rank RBPs by TPM (max TPM >= 10). """

def get_tpm(name,tpm):
    # Compute TPM metrics for RBP
    curr_rbp = tpm[tpm['gene_name']==name]
    id = curr_rbp.gene_id.values.tolist()
    result = []

    for i in id: # consider cases where multiple ids exist for RBP
        values = curr_rbp[curr_rbp['gene_id']==i].drop(['gene_id','gene_name'],axis=1).reset_index(drop=True).to_numpy()
        avg,min,max = 0,0,0
        if len(values) != 0:
            avg = np.mean(values)
            min = np.amin(values)
            max = np.amax(values)
        result.append([i,avg,min,max])

    if len(result) > 1: # more than 1 ID
        # compare average and return values with greater average
        high = 0
        index = 0
        for j in range(len(result)):
            if result[j][1] > high:
                high = result[j][1]
                index = j
        [id,avg, min, max] = result[index]
    else: # only one ID 
        [id,avg, min, max] = result[0]
    return id,avg, min, max

def filter_tpm(df):
    # Filter out RBPs by TPM
    # IF max TPM >= 10, keep RBP in list
    return df[df['max_TPM']>=10]

def rank_rbp(df):
    # Rank RBPs by TPM
    df['TPM_rank'] = df[['mean_TPM','max_TPM','min_TPM']].apply(tuple,axis=1).rank(method='dense',ascending=False).astype(int)
    df = df.sort_values('TPM_rank').reset_index(drop=True)
    return df

def main():
    rbp_file = snakemake.input[0] # RBP list (min 1 column: 'name')
    tpm_path = snakemake.params[0] # file path for coco_tpm.tsv
    outfile = snakemake.output[0] # file path for ranked RBPs

    # get TPM for cell lines not in tissue
    df_tpm = pd.read_csv(tpm_path,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    
    df_rbp = pd.read_csv(rbp_file,sep='\t',usecols=['name','protein_id'])
    
    df_exp_level = pd.DataFrame(columns=['id','name','mean_TPM','max_TPM','min_TPM'])

    for i in range(len(df_rbp)):
        curr_rbp = df_rbp.loc[i,'name']
        if curr_rbp == "TROVE2": # TROVE2 exists as RO60 in annotation (synonyms)
            curr_rbp = "RO60"
        id,avg, min, max = get_tpm(curr_rbp,df_tpm)
        result = pd.DataFrame({'id': id, 'name':curr_rbp,'protein_id':df_rbp.loc[i,'protein_id'],'mean_TPM':[avg],'max_TPM':[max],'min_TPM':[min]})
        df_exp_level = pd.concat([df_exp_level,result],ignore_index=True)

    # replace RO60 back to TROVE2
    df_exp_level = df_exp_level.replace("RO60","TROVE2")

    # Get RBP TPM info --> filter by TPM --> rank by TPM
    df_final = rank_rbp(filter_tpm(df_exp_level))
    df_final['type'] = 'RBP'
    
    # change column orders
    df_final = df_final[['id','name','type','protein_id','mean_TPM','max_TPM','min_TPM','TPM_rank']]

    df_final.to_csv(outfile,sep='\t',index=None)

if __name__ == '__main__':
    main()