rule snoglobe_uniq:
    """ Format snoglobe predictions for snoRNA-RBP p-val calculations. """
    input:
        os.path.join(config["data"]["snoglobe"],"pred_{sno}.98_3.gene.tsv")
    output:
        os.path.join(config["data"]["snoglobe_formatted"],"{sno}_uniq_regions.bed")
    params:
        config["data"]["snoglobe_formatted"]
    shell:
        "mkdir -p {params} && cut -f1-6 {input} | sort -k1,1 -k2,3n -u > {output}"

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

rule combine_snoglobe_htrri:
    message: "Combine snoGloBe predictions and HTRRI interactions."
    input:
        rules.snoglobe_uniq.output,
        rules.get_htrri.output
    output:
    shell: bedtools merge
