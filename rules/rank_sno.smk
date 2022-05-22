rule :
    input:
        config["data"]["snoRNA_list"]
    output:
    params:
        tpm_matrix = config["data"]["tpm"]
    shell:

"""
- Rank snoRNAs by expression level
- Then rank by # copies
"""