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