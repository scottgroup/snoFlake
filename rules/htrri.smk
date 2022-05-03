rule get_htrri:
    """ Extract snoRNA-snoRNA and snoRNA-RBP mRNA interactions from HTRRI dataset (PARIS, SPLASH, LIGR-seq). """
    input:
        htrri = config["data"]["HTRRI"],
        sno = config["data"]["snoRNA_list"],
        rbp = config["data"]["rbp_list"]
    output:
        os.path.join(config["outpath"],"filtered_HTRRI.tsv")
    shell:
        'python scripts/filter_htrri.py {input.htrri} {input.rbp} {input.sno} {output}'

