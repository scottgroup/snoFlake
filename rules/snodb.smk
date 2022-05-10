rule filter_snodb_host_gene:
    """ Find snoRNAs that have RBPs as host genes in snoDB """
    input:
        snodb = config["data"]["snoDB"],
        snolist = config["data"]["snoRNA_list"],
        rbplist = config["data"]["rbp_list"]
    output:
        os.path.join(config["outpath"],"snoDB_rbp_as_host_gene.tsv")
    shell:
        'python scripts/snoDB_host_gene.py {input.snodb} {input.snolist} {input.rbplist} {output}'