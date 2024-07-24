#!/usr/bin/env python

import pandas as pd
import itertools

"""
Find network motifs in snoFlake.
"""

def double_edges(df_edges):
    """
    Find all occurrences of double edges in the network: snoRNA & RBP target overlap AND RBP binds to snoRNA.
    """

    df_rbp_binds_to_sno = df_edges[df_edges['interaction']=='RBP_binds_to_snoRNA'].reset_index(drop=True)
    df_target_ovlp = df_edges[df_edges['interaction']=='sno_RBP_overlap'].reset_index(drop=True)

    df_target_ovlp['merged'] = df_target_ovlp['source'] + '+' + df_target_ovlp['target']
    df_rbp_binds_to_sno['merged'] = df_rbp_binds_to_sno['target'] + '+' + df_rbp_binds_to_sno['source']

    df_merged = pd.concat([df_target_ovlp,df_rbp_binds_to_sno],ignore_index=True)
    df_counts = df_merged['merged'].value_counts().reset_index()
    df_counts.columns = ['pairs','counts']
    df_counts[['snoRNA', 'RBP']] = df_counts['pairs'].str.split('+', n=1, expand=True)
    df_counts = df_counts[['snoRNA','RBP','counts']]

    df_double_edges = df_counts[df_counts['counts']==2].sort_values(by=['snoRNA']).reset_index(drop=True).drop(columns=['counts'])

    return df_double_edges


def search_motif(df,string):
    """
    Search for network motifs by grouping snoRNA-RBP interactions by snoRNA.
    """

    string_rbps = list(set(string['source'].to_list() + string['target'].to_list()))
    string['tuples'] = string[['source', 'target']].apply(tuple, axis=1)
    string_tuples_list = string['tuples'].values.tolist()
    df_by_sno = df.groupby('snoRNA', as_index=False).agg(list) # snoRNA [RBP list]

    two_rbps_count = 0

    df_out = pd.DataFrame(columns=['nodes'])

    for i in range(len(df_by_sno)):
        sno = df_by_sno.loc[i,'snoRNA']
        rbp_list = df_by_sno.loc[i,'RBP']

        # get list of RBPs that are part of STRING
        result = []
        for rbp in rbp_list:
            if rbp in string_rbps:
                result += [rbp]

        # check if RBP-RBP tuple present in STRING
        perm = itertools.permutations(result,2)
        string_tuples_result = []
        for p in list(perm):
            if p in string_tuples_list:
                string_tuples_result += [p]
                two_rbps_count += 1 # sno-RBP-RBP

        new_motif = [sno] + list(set(list(itertools.chain(*string_tuples_result))))

        if len(new_motif) > 1:
            df_out = pd.concat([df_out,pd.DataFrame([','.join(new_motif)],columns=df_out.columns)],ignore_index=True)

    # count motifs
    double_edge_count = len(df)
    print('Total number of 2-node graphlets: ' + str(double_edge_count))
    print('Total number of 3-node graphlets: ' + str(two_rbps_count))

    return df_out


def main():

    df_edges = pd.read_csv(snakemake.input.edges,sep='\t')
    
    db_edges = double_edges(df_edges)
    string = df_edges[df_edges['interaction']=='STRING'].reset_index(drop=True)
    string = string[['source','target']]

    search_motif(db_edges,string).to_csv(snakemake.output.network_motifs,sep='\t')


if __name__ == '__main__':
    main()