#!/usr/bin/env python

import pandas as pd
import py4cytoscape as p4c
import sys
import numpy as np


def node_weight(df):
    """
    Compute node size/weight based on TPM.
    """
    df['weight'] = 1
    df.loc[(df['max_TPM'] >= 10) & (df['max_TPM'] <100),'weight'] = 10
    df.loc[(df['max_TPM'] >= 100) & (df['max_TPM'] <1000),'weight'] = 100
    df.loc[df['max_TPM'] >= 1000,'weight'] = 1000
    return df


def edge_weight(df):
    """
    Compute edge thickness/weight based on p-values and scores.
    """
    df['weight'] = 1
    df.loc[(df['ENCODE_log_pval'] >= 10) & (df['ENCODE_log_pval'] < 20),'weight'] = 10
    df.loc[(df['ENCODE_log_pval'] >= 20) & (df['ENCODE_log_pval'] < 30),'weight'] = 100
    df.loc[df['ENCODE_log_pval'] >= 30,'weight'] = 1000
    df.loc[(df['sno_RBP_overlap_log_pval'] >= 10) & (df['sno_RBP_overlap_log_pval'] < 20),'weight'] = 10
    df.loc[(df['sno_RBP_overlap_log_pval'] >= 20) & (df['sno_RBP_overlap_log_pval'] < 30),'weight'] = 100
    df.loc[df['sno_RBP_overlap_log_pval'] >= 30,'weight'] = 1000
    df.loc[(df['STRING_score'] >= 925) & (df['STRING_score'] < 950),'weight'] = 10
    df.loc[(df['STRING_score'] >= 950) & (df['STRING_score'] < 975),'weight'] = 100
    df.loc[df['STRING_score'] >= 975,'weight'] = 1000

    # replace inf with 400 (sno_RBP_overlap)
    df.replace(np.inf, 400, inplace=True)

    return df


def search_network_motifs(motifs):
    """
    Extract double-edged network motifs using Cytoscape.
    """
    for i in range(len(motifs)):
        nodes_list = motifs.loc[i,'nodes'].split(',')
        motif_name = nodes_list[0] + '_' + str(len(nodes_list)) + '-node_graphlet'
        p4c.create_subnetwork(nodes=nodes_list, nodes_by_col='id', network="all_interactions", subnetwork_name=motif_name)
        p4c.clear_edge_bends(network=motif_name)
        p4c.layout_network('degree-circle', network=motif_name)

    return


def build_network(nodes,edges,motifs,smk_network):
    """
    Build snoFlake using Cytoscape --> Cytoscape app has to be OPEN
    """
    p4c.create_network_from_data_frames(nodes, edges, title="original", collection="snoFlake")

    # Create style
    style_name = "snoFlake Style"
    defaults = {'NODE_SHAPE':"ELLIPSE", 'NODE_TRANSPARENCY':255, 'NODE_BORDER_COLOR':"#252525", 'NODE_BORDER_WIDTH':3}
    nodeLabels = p4c.map_visual_property('node label', 'name', 'p')
    p4c.create_visual_style(style_name, defaults, [nodeLabels])
    p4c.set_visual_style(style_name)
    p4c.set_background_color_default('white', style_name='snoFlake Style')

    # Node design
    p4c.set_node_color_mapping('group', ['snoRNA', 'protein_coding'], ['#0C2C84', '#FECD87'], 'd', style_name='snoFlake Style')
    p4c.set_node_size_mapping('score', ['1','10','100','1000'], [40,60,80,100], 'd', style_name='snoFlake Style')
    p4c.set_node_label_color_mapping('group', ['snoRNA', 'protein_coding'], ['#0CE7D9', '#252525'], 'd', style_name='snoFlake Style')
    p4c.set_node_font_size_mapping('group', ['snoRNA', 'protein_coding'], [15, 12], 'd', style_name='snoFlake Style')

    # Edge design
    p4c.set_edge_color_mapping('interaction', ['RBP_binds_to_snoRNA', 'STRING', 'sno_RBP_overlap'], ['#01665E', '#FEC55F', '#9F9F9F'], 'd', style_name='snoFlake Style')
    p4c.set_edge_opacity_mapping('weight', ['1','10','100','1000'], [140,160,180,200], 'd', style_name='snoFlake Style')
    p4c.set_edge_line_width_mapping('weight', ['1','10','100','1000'], [2,4,6,8], 'd', style_name='snoFlake Style')

    # Network layout
    p4c.layout_network('degree-circle')

    # Filter out nodes with 0 edges
    p4c.create_degree_filter('degree filter', [1,5000], network="original")
    p4c.create_subnetwork(subnetwork_name='all_interactions')
    p4c.bundle_edges(network='all_interactions')
    p4c.clear_selection(network="original")

    # Find network motifs
    search_network_motifs(motifs)

    p4c.hide_panel('SOUTH')
    p4c.save_session(smk_network) # save session

    return


def main():
    
    # Snakemake variables
    smk_func = sys.argv[1]
    smk_nodes = sys.argv[2]
    smk_edges = sys.argv[3]


    if  smk_func == "data": # Get node and edge data

        smk_RBP_binds_to_sno = sys.argv[4]
        smk_STRING = sys.argv[5]
        smk_sno_RBP_overlap = sys.argv[6]
        smk_snoRNA = sys.argv[7]
        smk_RBP = sys.argv[8]

        # NODES
        sno = pd.read_csv(smk_snoRNA,sep='\t',usecols=['gene_id','gene_name','gene_biotype','min_TPM','max_TPM','rfam_id'])
        rbp = pd.read_csv(smk_RBP,sep='\t',usecols=['gene_name','gene_biotype','min_TPM','max_TPM'])
        rbp['gene_id'] = rbp['gene_name']
        rbp = rbp[['gene_id','gene_name','gene_biotype','min_TPM','max_TPM']]

        nodes = pd.concat([sno,rbp],ignore_index=True)
        nodes.gene_name.fillna(nodes.gene_id, inplace=True)
        nodes = node_weight(nodes)
        nodes.to_csv(smk_nodes,sep='\t',index=None)

        # EDGES
        rbp_binds_to_sno = pd.read_csv(smk_RBP_binds_to_sno,sep='\t')
        sno_rbp_ovlp = pd.read_csv(smk_sno_RBP_overlap,sep='\t')
        string = pd.read_csv(smk_STRING,sep='\t')
        edges = pd.concat([rbp_binds_to_sno,sno_rbp_ovlp,string],axis=0,ignore_index=True)
        edges = edges[['source','target','interaction','ENCODE_log_pval','sno_RBP_overlap_log_pval','STRING_score']]
        edges = edge_weight(edges)
        edges.to_csv(smk_edges,sep='\t',index=None)
    

    elif smk_func == "build": # Build snoFlake using Cytoscape --> Cytoscape app has to be OPEN

        smk_motifs = sys.argv[4]
        smk_network = sys.argv[5]

        nodes = pd.read_csv(smk_nodes,sep='\t')
        nodes.rename(columns={"gene_id":"id","gene_name":"name","gene_biotype":"group","weight":"score"},inplace=True)
        edges = pd.read_csv(smk_edges,sep='\t')
        motifs = pd.read_csv(smk_motifs,sep='\t')

        build_network(nodes,edges,motifs,smk_network)


if __name__ == '__main__':
    main()