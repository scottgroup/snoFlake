#!/usr/bin/env python

import py4cytoscape as p4c
import pandas as pd
import sys


def node_weights(df_gene):
    """
    Add node weights for size (based on max TPM).
    """
    df_gene['weight'] = 0 # arbitary value

    for i in range(len(df_gene)):
        max_tpm = df_gene.at[i,'max_TPM']
        if max_tpm < 1:
            df_gene.at[i,'weight'] = 0
        elif 1 <= max_tpm < 10:
            df_gene.at[i,'weight'] = 10
        elif 10 <= max_tpm < 100:
            df_gene.at[i,'weight'] = 100
        elif 100 <= max_tpm < 1000:
            df_gene.at[i,'weight'] = 1000
        elif 1000 <= max_tpm < 10000:
            df_gene.at[i,'weight'] = 10000
        else:
            df_gene.at[i,'weight'] = 100000
    
    return df_gene

"""
def default_settings():
    p4c.set_node_shape_default('ELLIPSE')

    column = 'type'
    values = ['snoRNA','RBP']
    colors = ['#FC9272','#A6BDDB']
    p4c.set_node_color_mapping(column,values,colors,mapping_type = 'discrete')
    # edge colors
    # STRING: '969696', rbp_sno_transcript: 'A6BDDB', sno_rbp_transcript: 'FDCDAC',sno_rbp_target_ovlp: '016450'
    return
"""

def main():

    # NODES
    sno = pd.read_csv(sys.argv[1],sep='\t',usecols=['gene_id','gene_name','gene_biotype','min_TPM','max_TPM'])
    rbp = pd.read_csv(sys.argv[2],sep='\t',usecols=['gene_name','gene_biotype','min_TPM','max_TPM'])
    rbp['gene_id'] = rbp['gene_name']
    rbp = rbp[['gene_id','gene_name','gene_biotype','min_TPM','max_TPM']]
    #pc = pd.read_csv(sys.argv[3],sep='\t',usecols=['gene_id','gene_name','gene_biotype'])
    nodes = pd.concat([node_weights(sno),node_weights(rbp)],ignore_index=True)
    nodes.gene_name.fillna(nodes.gene_id, inplace=True)
    nodes.to_csv(sys.argv[7],sep='\t',index=None)

    # EDGES
    rbp_binds_to_sno = pd.read_csv(sys.argv[6],sep='\t')
    sno_rbp_ovlp = pd.read_csv(sys.argv[5],sep='\t')
    string = pd.read_csv(sys.argv[4],sep='\t')
    edges = pd.concat([rbp_binds_to_sno,sno_rbp_ovlp,string],axis=0,ignore_index=True)
    edges.to_csv(sys.argv[8],sep='\t',index=None)

    # NETWORKS
    nodes_id = nodes.rename(columns={'gene_id':'id'})
    p4c.create_network_from_data_frames(nodes_id, edges, title="all_interactions", collection="sno_RBP_interaction_network")
    #default_settings()

    p4c.save_session(sys.argv[9])

if __name__ == '__main__':
    main()