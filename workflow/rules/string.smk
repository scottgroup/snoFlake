rule unzip:
    message: "Unzip STRING interactions .gz file."
    input:
        config["data"]["STRING"]
    output:
        temp(os.path.join(config["outpath"],"STRING_unzip.tsv"))
    shell:
        "gunzip -c {input} > {output}"

rule get_rbp_rbp:
    message: "Keep RBPs of interest only."
    input:
        string = rules.unzip.output,
        rbp = config["nodes"]["rbp_list_with_protein_id"]
    output:
        temp(os.path.join(config["temp"],"STRING_rbp_filtered.tsv"))
    shell:
        "python3 scripts/filter_STRING.py {input.rbp} {input.string} {output}"

rule get_rbp_rbp_count:
    message: "Count number of interactions left after filtering for RBP-RBP pairs."
    input:
        rules.get_rbp_rbp.output
    output:
        temp(os.path.join(config["temp"],"STRING_rbp_filtered_count.tsv"))
    shell:
        "int_cnt=$(grep -c ^ {input}) && int_cnt=$(($int_cnt-1)) && echo -e \"{rule}\t$int_cnt\" >> {output}"

rule filter_score:
    """ Keep STRING interactions above score threshold specified in config.json """
    input:
        string = rules.get_rbp_rbp.output
    output:
        temp(os.path.join(config["temp"],"STRING_score_filtered.tsv"))
    params:
        config["STRING_score_threshold"]
    shell: 
        "awk -F \'\t\' \'$4 > {params} {{print}}\' {input.string} > {output}"

rule filter_score_count:
    """ Count number of interactions left after filtering for interaction score """
    input:
        rules.filter_score.output
    output:
        temp(os.path.join(config["temp"],"STRING_score_filtered_count.tsv"))
    shell:
        "int_cnt=$(grep -c ^ {input}) && int_cnt=$(($int_cnt-1)) && echo -e \"{rule}\t$int_cnt\" >> {output}"

rule filter_physical_binding:
    """ Keep physical binding STRING interactions only """
    input:
        string = rules.filter_score.output
    output:
        os.path.join(config["outpath"],"STRING_physical_binding.tsv")
    shell:
        "echo -e \"RBP1\tRBP2\tinteraction\" >> {output} && "
        "awk -F \'\t\' \'$3 == \"binding\" {{print}}\' {input.string} | cut -f 1,2,5 | uniq >> {output}"

rule filter_physical_binding_count:
    """ Count number of interactions left after filtering for physical binding interactions """
    input:
        rules.filter_physical_binding.output
    output:
        temp(os.path.join(config["temp"],"STRING_physical_binding.tsv"))
    shell:
        "int_cnt=$(grep -c ^ {input}) && int_cnt=$(($int_cnt-1)) && echo -e \"{rule}\t$int_cnt\" >> {output}"

rule merge_all_STRING_counts:
    input:
        rbp = rules.get_rbp_rbp_count.output,
        score = rules.filter_score_count.output,
        phys = rules.filter_physical_binding_count.output
    output:
        temp(os.path.join(config["temp"],"STRING_count.tsv"))
    shell:
        "cat {input.rbp} > {output} && cat {input.score} >> {output} && cat {input.phys} >> {output}"