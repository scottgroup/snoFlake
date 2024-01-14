#!/usr/bin/python

import pandas as pd
import numpy as np
import math
pd.options.mode.chained_assignment = None  # default='warn'


""" Obtain snoRNA characteristics and create a ranking system. """
""" Box C/D count --> TPM filter --> # copies """


def sno_in_annotation(file):
    # extract snoRNAs from annotation
    df = pd.read_csv(file,sep='\t',usecols=['feature','gene_biotype','gene_name','gene_id'])
    df = df[df['feature']=='transcript']
    df = df[df['gene_biotype']=='snoRNA']
    df.drop(['feature', 'gene_biotype'], axis=1, inplace=True)
    df.drop_duplicates(inplace=True,ignore_index=True)
    return df

def filter_box(df,df_snodb):

    # filter out snoRNAs with weird gene_ids
    df = df[~df['id'].str.contains('cluster',na=False)].reset_index(drop=True)
    df = df[~df['id'].str.contains('snoDB',na=False)].reset_index(drop=True)

    # remove box H/ACA snoRNAs and scaRNAs
    for n in ['SNORA','SCA']:
        df = df[~df['name'].str.contains(n,na=False)].reset_index(drop=True)

    # check box type for snoRNAs with names not starting with 'SNORD'
    cd = df[df['name'].str.contains('SNORD',na=False)].reset_index(drop=True)
    others = df[~df['name'].str.contains('SNORD',na=False)].reset_index(drop=True)
    # check box type from snoDB for others
    for i in range(len(others)):
        id = others.loc[i,'id']
        info = df_snodb[df_snodb['ensembl_id'].str.contains(id,na=False)].reset_index(drop=True)
        if len(info) == 0:
            info = df_snodb[df_snodb['refseq_id']==id]
        if len(info) != 0 and info.loc[0,'box_type'] == 'C/D':
            curr_sno = pd.DataFrame({'id':[id],'name':[others.loc[i,'name']]})
            cd = pd.concat([cd,curr_sno],ignore_index=True)
    return cd

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

def filter_tpm(df,tpm):
    # filter snoRNAs found in annotation by TPM
    # get TPM for cell lines not in tissue
    df_tpm = pd.read_csv(tpm,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    # IF snoRNA TPM >= 10 in at least one cell line, THEN the snoRNA is kept in the list
    sno = pd.DataFrame(columns=['id','name'])
    for i in range(len(df)):
        _,_,max = get_tpm(df.loc[i,'id'],df_tpm)
        if max >= 10:
            curr_sno = pd.DataFrame({'id':[df.loc[i,'id']],'name':[df.loc[i,'name']]})
            sno = pd.concat([sno,curr_sno],ignore_index=True)
    return sno

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
    result = pd.DataFrame({'name': name, 'id':id,'snoDB_ID': snodb_id, 'rfam_id':rfam_id,'num_fam_members':[num_members],
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
    annotation = snakemake.input[0] # annotation tsv file  # snoRNA list (min 2 columns: 'id', 'name', ...)
    tpm_path = snakemake.params[0] # file path for coco_tpm.tsv
    snodb_path = snakemake.params[1] # file path for snoDB
    outfile = snakemake.output[0] # file path for ranked snoRNAs

    initial_sno = sno_in_annotation(annotation) # 1541 snoRNAs for hg38_Ensembl_V101_Scottlab_2020.tsv
    initial_sno = initial_sno.rename(columns={"gene_name": "name", "gene_id": "id"})

    # get TPM for cell lines not in tissue
    df_tpm = pd.read_csv(tpm_path,sep='\t',usecols=['gene_id','gene_name','HCT116_1','HCT116_2','MCF7_1','MCF7_2',
                                                    'PC3_1','PC3_2','SKOV_frg_1','SKOV_frg_2','TOV112D_1','TOV112D_2'])
    # snoDB in df
    df_snodb = pd.read_csv(snodb_path,sep='\t',usecols=['snodb_id','ensembl_id','rfam_id','refseq_id','gene_name','box_type','length','sequence'])

    box_cd = filter_box(initial_sno,df_snodb) # 566 snoRNAs
    tpm_10 = filter_tpm(box_cd,tpm_path) # 185 snoRNAs

    # extract 'id' and 'name' columns
    id_list = tpm_10['id'].values.tolist()
    name_list = tpm_10['name'].values.tolist()

    sno_char = pd.DataFrame(columns=['name','id','snoDB_ID','rfam_id','num_fam_members','num_highly_exp_fam_members','mean_TPM','max_TPM','min_TPM'])
    # get all characteristics for each snoRNA in the input list
    for i in range(len(id_list)):
        sno_char = pd.concat([sno_char,get_sno_info(id_list[i],name_list[i],df_tpm,df_snodb)],ignore_index=True)
    
    # Filter out by Rfam family member count (<10)
    rfam_filter = sno_char[sno_char['num_fam_members']<10] # 155 snoRNAs

    # rank snoRNAs
    df_final = rank_snornas(rfam_filter)

    # reformat df columns
    df_final['type'] = 'snoRNA'
    df_final = df_final[['id','name','type','snoDB_ID','rfam_id','num_fam_members','num_highly_exp_fam_members','mean_TPM','max_TPM','min_TPM']]    
    df_final['rank'] = df_final.reset_index().index + 1
    df_final.to_csv(outfile,sep='\t',index=None)

if __name__ == '__main__':
    main()