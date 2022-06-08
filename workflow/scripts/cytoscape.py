#!/usr/bin/env python

import py4cytoscape as p4c
import pandas as pd

def main():

    # NODES
    sno = pd.read_csv('../../results/snoRNA.tsv',sep='\t',usecols=['id','name','type'])
    rbp = pd.read_csv('../../results/rbp.tsv',sep='\t',usecols=['name','type'])
    # add additional 'id' column for RBP but with name
    rbp['id'] = rbp['name']
    rbp = rbp[['id','name','type']] # reorder df
    nodes = pd.concat([sno,rbp],ignore_index=True)

    # EDGES
    rbp_bind_to_sno = pd.read_csv('../../results/rbp_bind_to_sno.tsv',sep='\t')

    p4c.create_network_from_data_frames(nodes, rbp_bind_to_sno, title="my first network", collection="DataFrame Example")

if __name__ == '__main__':
    main()