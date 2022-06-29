#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def get_tpm(id,tpm):
    # Compute TPM metrics for snoRNA
    tpm_genes = tpm.iloc[:,0]

    # remove . for NR_
    id_format = id.split(".")
    id = id_format[0]

    values = tpm.loc[tpm_genes.isin([id]),].drop(['gene_id','gene_name'],axis=1).reset_index(drop=True).to_numpy()
    avg,min,max = 0,0,0
    if len(values) != 0:
        avg = np.mean(values)
        min = np.amin(values)
        max = np.amax(values)
    return avg, min, max

def main():
    df_tpm = pd.read_csv('../../resources/coco_tpm.tsv',sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    print(get_tpm('ENSG00000277194',df_tpm))

if __name__ == '__main__':
    main()