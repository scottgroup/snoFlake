rule unzip:
    """ Unzip STRING interactions .gz file """
    input:
        config["data"]["STRING"]
    output:
        temp(os.path.join(config["temp"],"STRING_unzip.tsv"))
    shell:
        "gunzip -c {input} > {output}"

rule get_rbp_rbp:
    """ Keep RBPs of interest only """
    input:
        string = rules.unzip.output,
        rbp = config["data"]["rbp_list_with_protein_id"]
    output:
        temp(os.path.join(config["temp"],"STRING_rbp_filtered.tsv"))
    shell:
        "python3 scripts/filter_STRING.py {input.rbp} {input.string} {output}"

rule filter_score:
    """ Keep STRING interactions above score threshold specified in config.json """
    input:
        rules.get_rbp_rbp.output
    output:
        temp(os.path.join(config["temp"],"STRING_score_filtered.tsv"))
    params:
        config["STRING_score_threshold"]
    shell: 
        "awk -F \'\t\' \'$7 > {params} {{print}}\' {input} > {output}"

rule filter_physical_binding:
    """ Keep physical binding STRING interactions only """
    input:
        rules.filter_score.output
    output:
        os.path.join(config["outpath"],"filtered_STRING.tsv")
    shell:
        "awk -F \'\t\' \'$3 == \"binding\" {{print}}\' {input} | cut -f 1,2 | uniq > {output}"