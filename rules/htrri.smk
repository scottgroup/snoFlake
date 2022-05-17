rule get_htrri:
    """ Extract snoRNA-snoRNA and snoRNA-RBP mRNA interactions from HTRRI dataset (PARIS, SPLASH, LIGR-seq). """
    input:
        htrri = config["data"]["HTRRI"],
        sno = config["data"]["snoRNA_list"],
        rbp = config["data"]["rbp_list"]
    output:
        interactions = os.path.join(config["outpath"],"filtered_HTRRI.tsv"),
        counts = temp(os.path.join(config["temp"],"filtered_HTRRI_count.tsv"))
    shell:
        'python scripts/filter_htrri.py {input.htrri} {input.rbp} {input.sno} {output} && '
        'int_cnt=$(grep -c ^ {output.interactions}) && int_cnt=$(($int_cnt-1)) && echo -e \"{rule}\t$int_cnt\" >> {output.counts}'