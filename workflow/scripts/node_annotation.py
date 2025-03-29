#!/usr/bin/env python

import pandas as pd

"""
Get snoRNA target info and RBP function and localization info.
"""

def reannotate_rbp(nodes,rbp_func,rbp_local):
    """
    Update rows with nan and all zeros.
    """
    nodes = nodes[nodes['gene_biotype']=='protein_coding']
    nodes = nodes[['gene_id']]

    # RBP function
    # Merge splicing-related columns
    rbp_func['Splicing'] = rbp_func[['Splicing regulation','Spliceosome']].max(axis=1)
    rbp_func.drop(columns=['Splicing regulation','Spliceosome'],inplace=True)

    # rename columns to remove space
    rbp_func = rbp_func.rename(columns={'RNA modification':'RNA_modification', "3' end processing":'3_end_processing', "rRNA processing":'rRNA_processing',
                                        'Ribosome & basic translation':'Ribosome_and_translation', 'RNA stability & decay':'RNA_stability_and_decay', 'microRNA processing':'microRNA_processing',
                                        'RNA localization':'RNA_localization', 'RNA export':'RNA_export', 'Translation regulation':'Translation_regulation', 'tRNA regulation':'tRNA_regulation',
                                        'mitochondrial RNA regulation':'mtRNA_regulation', 'Viral RNA regulation':'ViralRNA_regulation', 'snoRNA / snRNA / telomerase':'snoRNA_snRNA_telomerase',
                                        'P-body / stress granules':'P-body_stress_granules', 'Exon Junction Complex':'EJC'})
    cols = list(rbp_func.columns)
    cols.insert(1, cols.pop(-1))
    rbp_func = rbp_func[cols]

    rbp_func.loc[rbp_func['RBP'] == 'AARS', 'RBP'] = 'AARS1'
    rbp_func.loc[rbp_func['RBP'] == 'GARS', 'RBP'] = 'GARS1'
    rbp_func.loc[rbp_func['RBP'] == 'TROVE2', 'RBP'] = 'RO60'

    # add manual annotations (GO Biological Proccess) for RBPs that are missing annotations from ENCORE
    rbp_func.loc[rbp_func['RBP']=='APEX1'] = ['APEX1',0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='DDX43'] = ['DDX43',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='EIF4E'] = ['EIF4E',0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1]
    rbp_func.loc[rbp_func['RBP']=='ELAC2'] = ['ELAC2',0,0,1,0,0,1,0,0,0,0,1,1,0,0,0,0,0]
    rbp_func.loc[rbp_func['RBP']=='ELAVL1'] = ['ELAVL1',0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,1]
    rbp_func.loc[rbp_func['RBP']=='EXOSC10'] = ['EXOSC10',1,0,1,1,1,0,0,1,0,0,0,0,0,1,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='GARS1'] = ['GARS1',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='METTL1'] = ['METTL1',0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='MORC2'] = ['MORC2',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='RNF187'] = ['RNF187',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='RPS6'] = ['RPS6',0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='RYBP'] = ['RYBP',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='AATF'] = ['AATF',0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='FUBP3'] = ['FUBP3',0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='GRWD1'] = ['GRWD1',0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='SDAD1'] = ['SDAD1',0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='TBRG4'] = ['TBRG4',0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1]
    rbp_func.loc[rbp_func['RBP']=='UTP3'] = ['UTP3',0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1]

    mask = (rbp_func.iloc[:, 1:] == 0).all(axis=1)
    rbp_func.iloc[mask, -1] = 1

    # normalize RBP function matrix
    rbp_func.iloc[:, 1:] = rbp_func.iloc[:, 1:].div(rbp_func.iloc[:, 1:].sum(axis=1), axis=0)

    merged_func = pd.merge(nodes, rbp_func, left_on='gene_id', right_on='RBP', how='left')
    merged_func.drop(columns=['RBP'],inplace=True)

    # RBP localization
    merged_loc = pd.merge(nodes, rbp_local, left_on='gene_id', right_on='RBP', how='left')
    merged_loc.drop(columns=['RBP'],inplace=True)
    merged_loc = merged_loc.rename(columns={'PML bodies':'PML_bodies','Cajal bodies':'Cajal_bodies','P bodies':'P_bodies',
                                        'Nuclear release mitosis':'Nuclear_release_mitosis','Cell Cortex':'Cell_cortex'})

    # add manual annotations (GO subcellular location)
    merged_loc.loc[merged_loc['gene_id']=='AARS1'] = ['AARS1',1,0,0,0,0,1,1,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='ADAT1'] = ['ADAT1',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='AKAP8L'] = ['AKAP8L',1,0,1,1,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='APEX1'] = ['APEX1',1,1,1,0,0,1,1,0,0,1,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='APOBEC3C'] = ['APOBEC3C',1,0,0,0,0,1,0,0,1,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='DDX43'] = ['DDX43',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='DDX47'] = ['DDX47',1,1,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='DDX51'] = ['DDX51',1,1,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='EIF4E'] = ['EIF4E',1,0,1,0,0,1,0,0,1,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='ELAC2'] = ['ELAC2',1,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='ELAVL1'] = ['ELAVL1',1,0,0,0,0,1,0,0,1,1,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='EXOSC10'] = ['EXOSC10',1,1,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='EXOSC5'] = ['EXOSC5',1,1,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='FMR1'] = ['FMR1',1,1,0,0,1,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='FTO'] = ['FTO',1,0,1,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='GARS1'] = ['GARS1',0,0,0,0,0,1,1,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='GTF2F1'] = ['GTF2F1',1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='HNRNPA1'] = ['HNRNPA1',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='HNRNPL'] = ['HNRNPL',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='HNRNPM'] = ['HNRNPM',1,1,1,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='MBNL1'] = ['MBNL1',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='METTL1'] = ['METTL1',1,1,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='MORC2'] = ['MORC2',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='NIPBL'] = ['NIPBL',1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='POLR2G'] = ['POLR2G',1,0,0,0,0,0,0,0,1,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='RNF187'] = ['RNF187',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='RPS10'] = ['RPS10',0,1,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='RPS6'] = ['RPS6',1,1,0,0,0,1,0,0,0,1,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='RYBP'] = ['RYBP',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='SAFB'] = ['SAFB',1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='SDAD1'] = ['SDAD1',0,1,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='SF3B1'] = ['SF3B1',1,1,1,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='STAU2'] = ['STAU2',0,1,0,0,0,1,0,0,0,1,0,1,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='TIA1'] = ['TIA1',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='RO60'] = ['RO60',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='YWHAG'] = ['YWHAG',1,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    merged_loc.loc[merged_loc['gene_id']=='ZNF800'] = ['ZNF800',1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    # normalize RBP localization matrix
    merged_loc.iloc[:, 1:] = merged_loc.iloc[:, 1:].div(merged_loc.iloc[:, 1:].sum(axis=1), axis=0)
    merged_loc.iloc[:, 1:] = merged_loc.iloc[:, 1:].apply(lambda x: x.fillna(0) if x.isnull().all() else x, axis=1)

    return pd.merge(merged_func, merged_loc, on='gene_id', how='left')


def sno_target(nodes,rRNA,snRNA,snodb):
    """
    Check whether a snoRNA has a known canonical target or not. (1=present; 0=absent)
    """
    nodes = nodes[nodes['gene_biotype']=='snoRNA']
    nodes = nodes[['gene_id']]

    # snRNA targets
    targets = pd.merge(nodes, snRNA, left_on='gene_id', right_on='sno_id', how='left')
    targets = targets.rename(columns={'target': 'snRNA_target'})
    targets['snRNA_target'] = targets['snRNA_target'].apply(lambda x: 0 if pd.isna(x) else 1)
    targets.drop(columns=['sno_id'],inplace=True)

    # rRNA targets
    id = pd.merge(rRNA, snodb, left_on='snoDB_id', right_on='snodb_id', how='left')
    id = id.rename(columns={'ensembl_id': 'sno_id'})
    id['sno_id'] = id['sno_id'].fillna(id['snodb_id'])
    id.drop_duplicates(subset=['sno_id'],inplace=True)
    one_id = id[~id['sno_id'].str.contains(';', na=False)]
    two_ids = id[id['sno_id'].str.contains(';', na=False)]
    two_ids['sno_id'] = two_ids['sno_id'].str.split(';')
    two_ids[['A', 'B']] = pd.DataFrame(two_ids['sno_id'].tolist(), index=two_ids.index)
    two_ids1 = two_ids.copy()
    two_ids1['pair'] = two_ids1['A'] + '-' + two_ids1['B']
    two_ids2 = two_ids.copy()
    two_ids2['pair'] = two_ids2['B'] + '-' + two_ids2['A']
    two_ids = pd.concat([two_ids1, two_ids2], ignore_index=True)
    two_ids.drop(columns=['A','B','sno_id'],inplace=True)
    two_ids.rename(columns={'pair': 'sno_id'},inplace=True)
    id = pd.concat([one_id, two_ids], ignore_index=True)
    targets = pd.merge(targets, id, left_on='gene_id', right_on='sno_id', how='left')
    targets.rename(columns={'Type': 'rRNA_target'},inplace=True)
    targets.drop(columns=['snoDB_id','rRNA','snodb_id','sno_id'],inplace=True)
    targets['rRNA_target'] = targets['rRNA_target'].apply(lambda x: 0 if pd.isna(x) else 1)

    # canonical targets in total
    targets['canonical_target'] = targets[['snRNA_target','rRNA_target']].max(axis=1)

    return targets


def main():

    nodes = pd.read_csv(snakemake.input.nodes,sep='\t')
    rbp_func = pd.read_csv(snakemake.input.encore_matrix,usecols=['RBP','Splicing regulation','Spliceosome','RNA modification',"3' end processing","rRNA processing",
                                                'Ribosome & basic translation','RNA stability & decay','microRNA processing','RNA localization','RNA export',
                                                'Translation regulation','tRNA regulation','mitochondrial RNA regulation','Viral RNA regulation',
                                                'snoRNA / snRNA / telomerase','P-body / stress granules','Exon Junction Complex','Other'])
    rbp_local = pd.read_csv(snakemake.input.encore_matrix,usecols=['RBP','Nuclei','Nucleolus','Speckles','PML bodies',
                                                'Cajal bodies','Cytoplasm','Mitochondria','Golgi','P bodies',
                                                'ER','Cytoskeleton','Microtubule','Actin','Nuclear release mitosis','Cell Cortex'])
    sno_rRNA_target = pd.read_csv(snakemake.input.rRNA_targets,sep='\t',usecols=['snoDB_id','Type','rRNA'])
    sno_snRNA_target = pd.read_csv(snakemake.input.snRNA_targets,sep='\t',usecols=['sno_id','target'])
    snoDB = pd.read_csv(snakemake.input.snodb,sep='\t',usecols=['snodb_id','ensembl_id'])

    # RBP function and localization info
    reannotate_rbp(nodes,rbp_func,rbp_local).to_csv(snakemake.output.RBP_anno,sep='\t',index=None)
    # snoRNA target info
    sno_target(nodes,sno_rRNA_target,sno_snRNA_target,snoDB).to_csv(snakemake.output.sno_anno,sep='\t',index=None)


if __name__ == '__main__':
    main()