import pandas as pd

""" Get all snoRNAs with Ensembl id that have protein coding host genes in snoDB """

# Load snoDB
df = pd.read_csv('~/scratch/pipeline/data/snoDB.tsv',sep='\t',usecols=['ensembl_id','host_gene_name','host_biotype'])

# Format snoDB
df_pc = df[df['host_biotype']=='protein_coding'].dropna().reset_index()
df_pc.to_csv('~/scratch/pipeline/data/snoDB_pc_host_gene.tsv',sep='\t',columns=df.columns,header=True,index=False)
