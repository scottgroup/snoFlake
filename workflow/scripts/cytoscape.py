#!/usr/bin/env python

import py4cytoscape as p4c
import pandas as pd

def default_settings():
    p4c.set_node_shape_default('ELLIPSE')

    column = 'type'
    values = ['snoRNA','RBP']
    colors = ['#FC9272','#A6BDDB']
    p4c.set_node_color_mapping(column,values,colors,mapping_type = 'discrete')
    # edge colors
    # STRING: '969696', rbp_sno_transcript: 'A6BDDB', sno_rbp_transcript: 'FDCDAC',sno_rbp_target_ovlp: '016450'
    return

def main():

    # NODES
    sno = pd.read_csv('../../results/snoRNA.tsv',sep='\t',usecols=['id','name','type'])
    rbp = pd.read_csv('../../results/rbp.tsv',sep='\t',usecols=['name','type'])
    # add additional 'id' column for RBP but with name
    rbp['id'] = rbp['name']
    rbp = rbp[['id','name','type']] # reorder df
    nodes = pd.concat([sno,rbp],ignore_index=True)

    # EDGES
    rbp_bind_to_sno = pd.read_csv('../../results/rbp_bind_to_sno_transcript.tsv',sep='\t')
    sno_bind_to_rbp = pd.read_csv('../../results/sno_bind_to_rbp_transcript.tsv',sep='\t')
    sno_rbp_target_ovlp = pd.read_csv('../../results/significant_sno_rbp_target_overlaps.tsv',sep='\t')
    string = pd.read_csv('../../results/STRING_900_physical_binding.tsv',sep='\t')

    # NETWORKS
    p4c.create_network_from_data_frames(nodes, rbp_bind_to_sno, title="rbp_bind_to_sno", collection="sno_RBP_network")
    p4c.create_network_from_data_frames(nodes, sno_bind_to_rbp, title="sno_bind_to_rbp", collection="sno_RBP_network")
    p4c.create_network_from_data_frames(nodes, sno_rbp_target_ovlp, title="sno_rbp_target_ovlp", collection="sno_RBP_network")
    p4c.create_network_from_data_frames(nodes, string, title="STRING", collection="sno_RBP_network")

    #p4c.merge_networks(['rbp_bind_to_sno','sno_rbp_target_ovlp','sno_bind_to_rbp','STRING'],title='merged_network')
    default_settings()

if __name__ == '__main__':
    main()