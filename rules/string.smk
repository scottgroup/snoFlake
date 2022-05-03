rule string:
    """ Keep physical binding STRING interactions only """
    input:
        string = config["data"]["STRING"],
        rbp = config["data"]["rbp_list"]
    output:
        os.path.join(config["outpath"],"physical_binding_STRING.tsv")
    params:
        
    shell:
        "gunzip -c {input.string} | awk -F '\t' '($3 == "binding") && ($7 > 900) {print}' | cut -f 1,2 | uniq > $SCRATCH/p_3/STRING_temp.tsv"
