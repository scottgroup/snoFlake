rule rank_sno:
    message: "Obtain snoRNA characteristics and create a ranking system."
    input:
        config["data"]["annotation"]
    output:
        config["nodes"]["snoRNA_list"]
    params:
        tpm = config["data"]["tpm"],
        snodb = config["data"]["snoDB"]
    script:
        "../scripts/rank_sno.py"

rule rank_rbp:
    message: "Filter and rank RBPs by TPM."
    input:
        config["data"]["ENCODE_rbp_list"]
    output:
        config["nodes"]["rbp_list"]
    params:
        config["data"]["tpm"]
    script:
        "../scripts/rank_rbp.py"