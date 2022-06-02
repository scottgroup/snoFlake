#!/usr/bin/env python

import py4cytoscape as p4c
import pandas as pd

def main():
    #"https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/Overview-of-py4cytoscape.ipynb"
    edges = pd.read_csv('../../../network_output/rbp_bind_to_sno_transcript.tsv',sep='\t')
    nodes1 = pd.read_csv('../../../network_data/snoRNA.tsv',sep='\t')
    nodes2 = pd.read_csv('../../../network_data/rbp.tsv',sep='\t')
    nodes = pd.concat([nodes1,nodes2],ignore_index=True)


    p4c.create_network_from_data_frames(nodes1, edges, title="my first network", collection="DataFrame Example")

if __name__ == '__main__':
    main()