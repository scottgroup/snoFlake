rule get_sno:
    message: "Create a list of snoRNAs from the annotation and keep highly expressed ones."
    input:
        config["data"]["annotation"]
    output:
        config["data"]["snoRNA_list"]
    params:
        tpm = config["data"]["tpm"],
        snodb = config["data"]["snoDB"]
    shell:
        "python3 scripts/get_sno.py {input} {params.tpm} {params.snodb} {output}"

rule rank_sno:
    message: "Obtain snoRNA characteristics and create a ranking system."
    input:
        rules.get_sno.output
    output:
        os.path.join(config["outpath"],"snoRNA_ranking.tsv")
    params:
        tpm = config["data"]["tpm"],
        snodb = config["data"]["snoDB"]
    shell:
        "python3 scripts/rank_sno.py {input} {params.tpm} {params.snodb} {output}"

rule rank_rbp:
    message: "Get RBP TPM info and rank by TPM."
    input:
        config["data"]["ENCODE_rbp_list"]
    output:
        os.path.join(config["outpath"],"rbp_ranking.tsv")
    params:
        config["data"]["tpm"]
    shell:
        "python3 scripts/rank_rbp.py {input} {params} {output}"

rule format_rbp:
    message: "Format RBP file for Cytoscape integration."
    input:
        rules.rank_rbp.output
    output:
        config["data"]["rbp_list"]
    shell:
        "echo -e \"id\tname\ttype\" >> {output} && "
        "awk '{{print $1\"\t\"$2\"\t\"$7}}' {input} | sed 1d >> {output} "