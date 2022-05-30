#!/usr/bin/python

from cmath import nan
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import sys
import math

""" Obtain snoRNA characteristics and create a ranking system. """


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

def get_rfam_members(curr_rfam,snodb):
    # split curr_rfam by ; to detect more than one Rfam ID
    members = []
    if not curr_rfam != curr_rfam: # skip nan
        ids = curr_rfam.split(";")
        df_members = pd.DataFrame(columns=snodb.columns)
        # Get all family members by RfamID
        for id in ids:
            new_members = snodb.loc[snodb['rfam_id'].str.contains(id,na=False)].reset_index(drop=True)
            df_members = pd.concat([df_members,new_members],ignore_index=True)
        df_members = df_members.drop_duplicates(ignore_index=True) # remove members counted twice
        members = df_members[['ensembl_id','refseq_id','rfam_id','gene_name']].values.tolist()
    return members

def verify_rfam_members(members,tpm):
    # Check if family members are highly expressed (max TPM >= 10)
    exp_members = []
    for mem in members:
        id = mem[0]
        if id != id: # check for nan in ensembl_id column
            id = mem[1]
        if id != id: # no ensembl_id nor refseq_id
            continue
        # check for ; in id
        indiv_id = id.split(";")
        for i in indiv_id:
            _,_,max = get_tpm(i,tpm)
            if max >= 10: ### CAN CHANGE TPM THRESHOLD HERE
                exp_members.append(i)
    return exp_members

def get_sno_info(id,name,tpm,snodb):
    # Compute TPM values (AVG, MIN, MAX) for specified snoRNA.
    avg, min, max = get_tpm(id,tpm)
    ### COMMENT OUT THIS SECTION IF EXACT VALUES NEEDED
    # round numbers 
    rounded = []
    for num in [avg,min,max]:
        if num > 0:
            dig = int(math.log10(num))
            num = np.round(num,decimals=-dig)
        rounded.append(num)
    avg = rounded[0]
    min = rounded[1]
    max = rounded[2]

    # Get other snoRNA characteristics including number of family members based on RfamID
    sno_info = snodb.loc[snodb['ensembl_id'].str.contains(id,na=False)].reset_index(drop=True)
    if sno_info.shape[0] > 0 :
        snodb_id = sno_info.at[0,'snodb_id']
        rfam_id = sno_info.at[0,'rfam_id']
    else: # snoRNAs with id starting with "NR_"
        snodb_id = snodb.loc[snodb['refseq_id'].str.contains(id,na=False)].reset_index(drop=True).at[0,'snodb_id']
        rfam_id = snodb.loc[snodb['refseq_id'].str.contains(id,na=False)].reset_index(drop=True).at[0,'rfam_id']
    fam_members_list = get_rfam_members(rfam_id,snodb)
    highly_exp_members = len(verify_rfam_members(fam_members_list,tpm))
    num_members = len(fam_members_list)

    # 0 --> no Rfam ID present for that snoRNA
    result = pd.DataFrame({'snoRNA': name, 'Ensembl_ID/RefSeq_ID':id,'snoDB_ID': snodb_id, 'rfam_id':rfam_id,'num_fam_members':[num_members],
                            'num_highly_exp_fam_members':[highly_exp_members],'mean_TPM':[avg],'max_TPM':[max],'min_TPM':[min]})
    return result

def rank_snornas(df):
    # group by number of highly expressed Rfam family members --> rank within group by TPM
    # get member count as list (for grouping)
    num_members_list = list(df.groupby('num_highly_exp_fam_members').groups.keys())

    # rank by TPM values for each num_members group
    dic_df = {}
    num_members_list.remove(0) # combine group 0 and 1
    for num in num_members_list:
        df_temp = df.loc[df['num_highly_exp_fam_members'] == num]
        if num == 1:
            df_temp = pd.concat([df_temp,df.loc[df['num_highly_exp_fam_members'] == 0]],ignore_index=True)
        df_temp['TPM_rank'] = df_temp[['mean_TPM','max_TPM','min_TPM']].apply(tuple,axis=1).rank(method='dense',ascending=False).astype(int)
        df_temp = df_temp.sort_values('TPM_rank').reset_index(drop=True)
        dic_df[str(num)] = df_temp

    # concat groups into one df
    df_result = pd.DataFrame(columns=dic_df['1'].columns)
    for num in num_members_list:
        df_result = pd.concat([df_result,dic_df[str(num)]],ignore_index=True)

    return df_result

def main():
    sno_file = sys.argv[1] # snoRNA list (min 2 columns: 'id', 'name', ...)
    tpm_path = sys.argv[2] # file path for coco_tpm.tsv
    snodb_path = sys.argv[3] # file path for snoDB
    outfile = sys.argv[4] # file path for ranked snoRNAs

    # get TPM for cell lines not in tissue
    df_tpm = pd.read_csv(tpm_path,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    # snoDB in df
    df_snodb = pd.read_csv(snodb_path,sep='\t',usecols=['snodb_id','ensembl_id','rfam_id','refseq_id','gene_name','box_type','length','sequence'])

    df_sno = pd.read_csv(sno_file,sep='\t',usecols=['id','name']) 
    # extract 'id' and 'name' columns
    id_list = df_sno['id'].values.tolist()
    name_list = df_sno['name'].values.tolist()

    sno_char = pd.DataFrame(columns=['snoRNA','Ensembl_ID/RefSeq_ID','snoDB_ID','rfam_id','num_fam_members','num_highly_exp_fam_members','mean_TPM','max_TPM','min_TPM'])
    # get all characteristics for each snoRNA in the input list
    for i in range(len(id_list)):
        sno_char = pd.concat([sno_char,get_sno_info(id_list[i],name_list[i],df_tpm,df_snodb)],ignore_index=True)
    
    # ONLY WANT SNO INFO
    #sno_char.to_csv(outfile,sep='\t',index=None)
    # RANK SNO
    rank_snornas(sno_char).to_csv(outfile,sep='\t',index=None)

if __name__ == '__main__':
    main()