rule get_rbp_rbp:
    """ Keep physical binding STRING interactions only """
    input:
        string = config["data"]["STRING"],
        rbp = config["data"]["rbp_list"]
    output:
        os.path.join(config["outpath"],"physical_binding_STRING.tsv")
    params:
        
    shell:
        "gunzip -c {input.string} | awk -F '\t' '($3 == "binding") && ($7 > 900) {print}' | cut -f 1,2 | uniq > $SCRATCH/p_3/STRING_temp.tsv"

rule filter_score:
    """ Keep STRING interactions above score threshold specified in config.json """
    input:
        rules.get_rbp_rbp.output
    params:
        config["STRING_score_threshold"]
    shell: 


rule filter_physical_binding:
    """ Keep physical binding STRING interactions only """
    input:
        rules.filter_score.output
    output:
        os.path.join(config["outpath"],"filtered_STRING.tsv")