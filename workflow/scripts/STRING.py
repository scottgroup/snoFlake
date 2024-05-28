#!/usr/bin/env python


""" Filter STRING interactions by only keeping RBP-RBP pairs specified in the RBP list """

import pandas as pd
from composite_score_ENCODE import min_max_normalize


def select_interactions(interactions_df,info_df,rbp_list_df):
    """
    Select RBP-RBP physical binding interactions of interest.
    """

    # Replace protein ids with names
    for i in range(len(rbp_list_df)):
        name = rbp_list_df.loc[i,'gene_name']
        id = info_df[info_df['preferred_name']==name].reset_index(drop=True).loc[0,'string_protein_id']
        interactions_df.replace(id,name,inplace=True)

    # Drop rows that contain protein ids
    interactions_df = interactions_df[interactions_df['protein1'].str.contains("9606.ENSP")==False]
    interactions_df = interactions_df[interactions_df['protein2'].str.contains("9606.ENSP")==False].reset_index(drop=True)

    return interactions_df


def remove_duplicate_edges(interactions_df):
    """
    Remove duplicate edges by ignoring edge direction.
    """
    
    # Store interaction pairs as sets (unordered)
    final_interactions_list = []
    interactions_list = interactions_df.values.tolist()
    final_interactions_df = pd.DataFrame(columns=interactions_df.columns)

    for row in interactions_list:
        curr_set = set(row)
        if curr_set not in final_interactions_list:
            final_interactions_list.append(curr_set)
            curr_df = pd.DataFrame([row],columns=final_interactions_df.columns)
            final_interactions_df = pd.concat([final_interactions_df,curr_df],ignore_index=True)

    final_interactions_df.rename(columns={"protein1": "source", "protein2": "target"},inplace=True)
    final_interactions_df['interaction'] = "STRING"
    return final_interactions_df


def main():
    interactions = snakemake.input.interactions # STRING interactions above the combined score threshold
    info = snakemake.input.info # List of proteins and their ids from the STRING database
    rbp_list = snakemake.params.RBP_list # List of RBPs to be included in the network
    outfile = snakemake.output[0]

    interactions_df = pd.read_csv(interactions,delim_whitespace=True,usecols=['protein1','protein2','combined_score'])
    interactions_df['protein1'] = interactions_df['protein1'].astype(str)
    interactions_df['protein2'] = interactions_df['protein2'].astype(str)
    info_df = pd.read_csv(info,sep='\t',usecols=['#string_protein_id','preferred_name'])
    info_df.rename(columns={"#string_protein_id": "string_protein_id"},inplace=True)
    rbp_list_df = pd.read_csv(rbp_list,sep='\t',usecols=['gene_name'])

    selected_interactions_df = select_interactions(interactions_df,info_df,rbp_list_df)
    selected_interactions_df = remove_duplicate_edges(selected_interactions_df)

    # normalize score
    selected_interactions_df['normalized_score'] = min_max_normalize(selected_interactions_df['combined_score'])
    selected_interactions_df = selected_interactions_df[['source','target','normalized_score','interaction']]
    selected_interactions_df.to_csv(outfile, sep='\t',index=False)



if __name__ == '__main__':
    main()