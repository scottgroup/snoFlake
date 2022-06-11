rule get_rbp_rbp:
    message: "Keep RBPs of interest only."
    input:
        config["data"]["STRING"],
        config["data"]["ENCODE_rbp_list"]
    output:
        temp(os.path.join(config["outpath"],"STRING_rbp_filtered.tsv"))
    script:
        "../scripts/filter_STRING.py"

rule filter_STRING:
    message: "Keep physical binding interactions that are above the threshold."
    input:
        rules.get_rbp_rbp.output
    output:
        os.path.join(config["outpath"],"STRING_900_physical_binding.tsv")
    params:
        config["filters"]["STRING_score_threshold"]
    shell: 
        "echo -e \"source\ttarget\tinteraction\" > {output} && "
        "awk -F \'\t\' \'($4 > {params} && $3 == \"binding\") {{print}}\' {input} | cut -f 1,2,5 | sort -u >> {output}"