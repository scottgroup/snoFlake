#!/usr/bin/env python

""" snoRNA and RBP overlapping target analysis."""

import sys
import os
import pandas as pd

def main():
    fpath = sys.argv[1] # path to directory containing target lists
    #pairs = sys.argv[2] # file containing interested snoRNA-RBP pairs
    outfile = sys.argv[2] # file containing output interactions for cytoscape

    #pairs_list = pd.read_csv(pairs,sep='\t',header=None).values.tolist()

    sno = 'ENSG00000277194'
    rbp = ['EFTUD2','PRPF8','SERBP1']

    # get target lists for each RBP
    d = {}
    for r in rbp:
        infile = os.path.join(fpath, '%s_%s_targetlist.txt' % (sno,r))
        open_file = open(infile, "r")
        target_list = open_file.read().split(",")
        d[r] = target_list
        print(r)
        print(len(d[r]))
        open_file.close()
    
    eftud2_prpf8 = set(d['EFTUD2']).intersection(d['PRPF8']) #34
    eftud2_serbp1 = set(d['EFTUD2']).intersection(d['SERBP1']) #0
    prpf8_serbp1 = set(d['PRPF8']).intersection(d['SERBP1']) #1
    
    for i in [eftud2_prpf8,eftud2_serbp1,prpf8_serbp1]:
        print(i)
    
    # create network file
    interaction = pd.DataFrame(columns=['source','target','interaction'])

    # EFTUD2 & PRPF8
    ep_sno = pd.DataFrame(columns=['source','target'])
    ep_sno['target'] = list(eftud2_prpf8)
    ep_sno['source'] = sno
    ep_eftud2 = pd.DataFrame(columns=['source','target'])
    ep_eftud2['target'] = list(eftud2_prpf8)
    ep_eftud2['source'] = 'EFTUD2'
    ep_prpf8 = pd.DataFrame(columns=['source','target'])
    ep_prpf8['target'] = list(eftud2_prpf8)
    ep_prpf8['source'] = 'PRPF8'
    ep = pd.concat([ep_sno,ep_eftud2,ep_prpf8],ignore_index=True)
    ep['interaction'] = 'EFTUD2_PRPF8'

    # PRPF8 & SERBP1
    data = [[sno, 'ENSG00000074800'], ['PRPF8', 'ENSG00000074800'], ['SERBP1', 'ENSG00000074800']]
    ps = pd.DataFrame(data, columns=['source', 'target'])
    ps['interaction'] = 'PRPF8_SERBP1'

    # targets that are unique to RBPs
    for a in list(eftud2_prpf8):
        d['EFTUD2'].remove(a)
        d['PRPF8'].remove(a)
    d['SERBP1'].remove('ENSG00000074800')
    d['PRPF8'].remove('ENSG00000074800')

    eftud2_sno = pd.DataFrame(columns=['source','target'])
    eftud2_sno['target'] = d['EFTUD2']
    eftud2_sno['source'] = sno
    eftud2_rbp = pd.DataFrame(columns=['source','target'])
    eftud2_rbp['target'] = d['EFTUD2']
    eftud2_rbp['source'] = 'EFTUD2'
    eftud2 = pd.concat([eftud2_sno,eftud2_rbp],ignore_index=True)
    eftud2['interaction'] = 'EFTUD2'

    prpf8_sno = pd.DataFrame(columns=['source','target'])
    prpf8_sno['target'] = d['PRPF8']
    prpf8_sno['source'] = sno
    prpf8_rbp = pd.DataFrame(columns=['source','target'])
    prpf8_rbp['target'] = d['PRPF8']
    prpf8_rbp['source'] = 'PRPF8'
    prpf8 = pd.concat([prpf8_sno,prpf8_rbp],ignore_index=True)
    prpf8['interaction'] = 'PRPF8'

    serbp1_sno = pd.DataFrame(columns=['source','target'])
    serbp1_sno['target'] = d['SERBP1']
    serbp1_sno['source'] = sno
    serbp1_rbp = pd.DataFrame(columns=['source','target'])
    serbp1_rbp['target'] = d['SERBP1']
    serbp1_rbp['source'] = 'SERBP1'
    serbp1 = pd.concat([serbp1_sno,serbp1_rbp],ignore_index=True)
    serbp1['interaction'] = 'SERBP1'

    interaction = pd.concat([ep,ps,eftud2,prpf8,serbp1],ignore_index=True)
    interaction.to_csv(outfile,sep='\t',index=None)

if __name__ == '__main__':
    main()