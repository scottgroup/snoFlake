#!/usr/bin/python

import pandas as pd
import numpy as np
import sys

def blastn(copy_list):
    # Compare % identity between snoRNA copies
    return

def get_copies(id,curr_rfam,snodb):
    # Get all snoRNA copies by RfamID
    df_copies = snodb.loc[snodb['rfam_id'] == curr_rfam]
    copies = df_copies[['ensembl_id','refseq_id','gene_name']].values.tolist()
    for i in copies[:-1]:
        if i[0] == id or i[1] == id:
            copies.remove(i)
    return copies

def get_sno_info(id,name,tpm,snodb):
    # Compute TPM values (AVG, MIN, MAX) for specified snoRNA.
    tpm_genes = tpm.iloc[:,0]
    values = tpm.loc[tpm_genes.isin([id]),].drop(['gene_id','gene_name'],axis=1).reset_index(drop=True).to_numpy()
    avg = np.mean(values)
    min = np.amin(values)
    max = np.amax(values)
    # Get other snoRNA characteristics including number of copies based on RfamID
    snodb_id = snodb.loc[snodb['ensembl_id'].str.contains(id,na=False)].reset_index(drop=True)
    rfam_id = snodb.loc[snodb['ensembl_id'].str.contains(id,na=False)].reset_index(drop=True)
    if snodb_id.shape[0] > 0 :
        snodb_id = snodb_id.at[0,'snodb_id']
        rfam_id = rfam_id.at[0,'rfam_id']
    else: # snoRNAs with id starting with "NR_"
        snodb_id = snodb.loc[snodb['refseq_id'].str.contains(id,na=False)].reset_index(drop=True).at[0,'snodb_id']
        rfam_id = snodb.loc[snodb['refseq_id'].str.contains(id,na=False)].reset_index(drop=True).at[0,'rfam_id']
    copies_list = get_copies(id,rfam_id,snodb)
    result = pd.DataFrame({'snoRNA': name, 'Ensembl_ID/RefSeq_ID':id,'snoDB_ID': snodb_id, 'num_copies':[len(copies_list)],'mean_TPM':[avg],'max_TPM':[max],'min_TPM':[min]})
    return result

def main():
    sno_file = sys.argv[1] # snoRNA list (min 2 columns: 'id', 'name', ...)
    tpm_path = sys.argv[2] # file path for coco_tpm.tsv
    snodb_path = sys.argv[3] # file path for snoDB

    # get TPM for cell lines not in tissue
    df_tpm = pd.read_csv(tpm_path,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    # snoDB in df
    df_snodb = pd.read_csv(snodb_path,sep='\t',usecols=['snodb_id','ensembl_id','rfam_id','refseq_id','gene_name','box_type','length','sequence'])

    df_sno = pd.read_csv(sno_file,sep='\t') 
    # extract 'id' and 'name' columns
    id_list = df_sno['id'].values.tolist()
    name_list = df_sno['name'].values.tolist()

    sno_char = pd.DataFrame(columns=['snoRNA','Ensembl_ID/RefSeq_ID','snoDB_ID','num_copies','mean_TPM','max_TPM','min_TPM'])
    # get all characteristics for each snoRNA in the input list
    for i in range(len(id_list)):
        sno_char = pd.concat([sno_char,get_sno_info(id_list[i],name_list[i],df_tpm,df_snodb)],ignore_index=True)

    # rank snoRNAs by their characteristics

if __name__ == '__main__':
    main()