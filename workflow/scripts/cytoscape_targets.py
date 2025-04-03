#!/usr/bin/env python

import pandas as pd
import os

"""
Reformat significant snoRNA-RBP overlapping targets for Cytoscape.
"""


def merge_sig_targets(path,edges):
    df_sig = pd.DataFrame(columns=['snoRNA','RBP','target_id','num_ovlp_nt'])
    for i in range(len(edges)):
        sno = edges.loc[i,'source']
        rbp = edges.loc[i,'target']
        ovlp = pd.read_csv(os.path.join(path,sno+'_RBP.tsv'),sep='\t')
        ovlp = ovlp[ovlp['RBP']==rbp]
        ovlp = ovlp[['snoRNA','RBP','target_id','num_ovlp_nt']]
        ovlp.sort_values(by=['snoRNA','RBP','target_id','num_ovlp_nt'],inplace=True,ascending=[False,False,False,False])
        ovlp.drop_duplicates(subset=['snoRNA','RBP','target_id'],keep="first",inplace=True) # remove duplicated targets and keep first occurence
        df_sig = pd.concat([df_sig,ovlp],ignore_index=True)
    return df_sig


def create_nodes(df,exp_pc_genes,nodes):
    # If target corresponds to an RBP in the network, it will not be included in this node file, but will be part of the edge file.
    df_nodes = pd.DataFrame()

    df.drop_duplicates(subset=['target_id'],inplace=True)
    df = df.reset_index(drop=True)
    df_nodes['gene_id'] = df['target_id'].values.tolist()
    df_nodes = df_nodes.merge(exp_pc_genes,on='gene_id',how='left')
    df_nodes = df_nodes.dropna(subset=["max_TPM"])# remove lowly expressed genes
    df_nodes = df_nodes.merge(nodes,left_on='gene_name',right_on='gene_id',how='left')
    df_nodes = df_nodes[df_nodes["weight"].isna()] # remove RBPs
    df_nodes['gene_biotype'] = 'protein_coding_target'
    df_nodes = df_nodes[['gene_id_x','gene_name_x','gene_biotype','avg_TPM','min_TPM_x','max_TPM_x']]
    df_nodes.rename(columns={'gene_id_x':'gene_id','gene_name_x':'gene_name',
                            'min_TPM_x':'min_TPM','max_TPM_x':'max_TPM'},inplace=True)
    df_nodes["gene_name"] = df_nodes["gene_name"].fillna(df_nodes["gene_id"])
    # weights
    df_nodes['weight'] = 1
    df_nodes.loc[(df_nodes['max_TPM'] >= 10) & (df_nodes['max_TPM'] < 100),'weight'] = 10
    df_nodes.loc[(df_nodes['max_TPM'] >= 100) & (df_nodes['max_TPM'] < 1000),'weight'] = 100
    df_nodes.loc[df_nodes['max_TPM'] >= 1000,'weight'] = 1000
    return df_nodes


def create_edges(df,exp_pc_genes,nodes):
    df = df.merge(exp_pc_genes,left_on='target_id',right_on='gene_id',how='left')
    df = df[['snoRNA','RBP','target_id','num_ovlp_nt','gene_name']]
    df = df.merge(nodes,on='gene_name',how='left')
    df.loc[df["gene_id"].notna(), "target_id"] = df["gene_id"]
    df = df[['snoRNA','RBP','target_id','num_ovlp_nt']]
    # weights
    df['weight'] = 0.1
    df.loc[(df['num_ovlp_nt'] >= 10) & (df['num_ovlp_nt'] < 20),'weight'] = 0.2
    df.loc[(df['num_ovlp_nt'] >= 20) & (df['num_ovlp_nt'] < 30),'weight'] = 0.3
    df.loc[df['num_ovlp_nt'] >= 30,'weight'] = 0.4
    # gather interactions
    df_edges_sno = pd.DataFrame(columns=['source','target','interaction','num_ovlp_nt','weight'])
    df_edges_sno['source'] = df['snoRNA']
    df_edges_sno['target'] = df['target_id']
    df_edges_sno['interaction'] = 'sno_RBP_ovlp_target'
    df_edges_sno['num_ovlp_nt'] = df['num_ovlp_nt']
    df_edges_sno['weight'] = df['weight']
    df_edges_rbp = pd.DataFrame(columns=['source','target','interaction','num_ovlp_nt','weight'])
    df_edges_rbp['source'] = df['RBP']
    df_edges_rbp['target'] = df['target_id']
    df_edges_rbp['interaction'] = 'sno_RBP_ovlp_target'
    df_edges_rbp['num_ovlp_nt'] = df['num_ovlp_nt']
    df_edges_rbp['weight'] = df['weight']
    df_edges = pd.concat([df_edges_sno,df_edges_rbp],ignore_index=True)
    return df_edges


def main():

    target_path = snakemake.params.target_path
    edges = pd.read_csv(snakemake.input.edges,sep='\t')
    edges = edges[edges['interaction']=='sno_RBP_overlap'].reset_index(drop=True)
    exp_pc_genes = pd.read_csv(snakemake.input.exp_pc_genes,sep='\t',usecols=['gene_id','gene_name','gene_biotype','avg_TPM','min_TPM','max_TPM'])
    nodes = pd.read_csv(snakemake.input.nodes,sep='\t')
    nodes = nodes[nodes['gene_biotype']=='protein_coding']

    sig_targets = merge_sig_targets(target_path,edges)
    create_nodes(sig_targets,exp_pc_genes,nodes).to_csv(snakemake.output.nodes,sep='\t',index=None)
    create_edges(sig_targets,exp_pc_genes,nodes).to_csv(snakemake.output.edges,sep='\t',index=None)
 

if __name__ == '__main__':
    main()