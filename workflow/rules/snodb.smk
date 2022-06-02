rule filter_snodb_host_gene:
    """ Find snoRNAs that have RBPs as host genes in snoDB """
    input:
        snodb = config["data"]["snoDB"],
        snolist = config["data"]["snoRNA_list"],
        rbplist = config["data"]["rbp_list"]
    output:
        interactions = os.path.join(config["outpath"],"snoDB_rbp_as_host_gene.tsv"),
        counts = temp(os.path.join(config["temp"],"snoDB_rbp_as_host_gene_count.tsv"))
    shell:
        'python scripts/snoDB_host_gene.py {input.snodb} {input.snolist} {input.rbplist} {output.interactions} && '
        'int_cnt=$(grep -c ^ {output.interactions}) && int_cnt=$(($int_cnt-1)) && echo -e \"{rule}\t$int_cnt\" >> {output.counts}'