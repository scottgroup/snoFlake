#!/usr/bin/python

import pandas as pd
import numpy as np
import sys

""" Compute TPM values (AVG, MIN, MAX) for specified snoRNA. """

def compute(id,tpm):
    tpm_genes = tpm.iloc[:,0]
    values = tpm.loc[tpm_genes.isin([id]),].drop(['gene_id','gene_name'],axis=1).reset_index(drop=True).to_numpy()
    avg = np.mean(values)
    min = np.amin(values)
    max = np.amax(values)
    return avg,min,max

def main():
    sno_list = sys.argv[1] # snoRNA list
    tpm_path = sys.argv[2] # file path for coco_tpm.tsv

    # get TPM for cell lines not in tissue
    tpm_df = pd.read_csv(tpm_path,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])

    # get mean, min and max TPM for snoRNA
    compute(gene_id,tpm_df)

if __name__ == '__main__':
    main()