rule access_snodb:
    """ Get all snoRNAs with Ensembl id that have protein coding host genes in snoDB """
    input:
        config["data"]["snoDB"]
    output:
        temp(os.path.join(config["temp"],"snodb_pc_host_gene.tsv"))
    shell:
        'python scripts/get_snoDB_pc_host_gene.py '

rule filter_host_gene:
    input:
        rules.access_snodb.output
    output: