#!/usr/bin/python

import pandas as pd
import numpy as np
import sys

def blastn(copy_list):
    # Compare % identity between snoRNA copies
    return

def get_copies(id,snodb):
    # Get all snoRNA copies by RfamID
    copies = []
    return

def get_sno_info(sno,tpm,snodb):
    
    # Compute TPM values (AVG, MIN, MAX) for specified snoRNA.
    tpm_genes = tpm.iloc[:,0]
    values = tpm.loc[tpm_genes.isin([id]),].drop(['gene_id','gene_name'],axis=1).reset_index(drop=True).to_numpy()
    avg = np.mean(values)
    min = np.amin(values)
    max = np.amax(values)
    # Get other snoRNA characteristics including number of copies based on RfamID
    result = pd.DataFrame({'snoRNA': 0, 'snoDB_ID': 2, 'Ensembl_ID':3,'num_copies':x,'mean_TPM':avg,'max_TPM':max,'min_TPM':min})
    return result

def main():
    sno_list = sys.argv[1] # snoRNA list 
    tpm_path = sys.argv[2] # file path for coco_tpm.tsv
    snodb_path = sys.argv[3] # file path for snoDB

    # get TPM for cell lines not in tissue
    df_tpm = pd.read_csv(tpm_path,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    # snoDB in df
    df_snodb = pd.read_csv(snodb_path,sep='\t',usecols=['snodb_id','ensembl_id','rfam_id','gene_name','box_type','length','sequence'])

    # get all characteristics for each snoRNA in the input list
    sno_list = pd.read_csv(sno_list,sep='\t',header=None).iloc[:,0].values.tolist()
    for sno in sno_list:
        get_sno_info(sno,df_tpm,df_snodb)

if __name__ == '__main__':
    main()