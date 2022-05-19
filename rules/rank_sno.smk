rule :
    input:
        config["data"]["snoRNA_list"]
    output:
    params:
        tpm_matrix = config["data"]["tpm"]
    shell:
