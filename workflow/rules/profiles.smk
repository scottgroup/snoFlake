rule get_fasta:
    message: "Create snoRNA fasta file for interaction_region.py."
    input:
    output:
    shell:

rule interaction_region:
    message: "Plot interaction region profile in snoRNAs."
    input:
        rules.get_fasta.output
    output:
    script:
        "../scripts/merge_snoglobe_htrri.py"

rule get_chromo_info:
    message: "Create chromo.txt file for nts_flanking_exons.py."


rule nts_flanking_exons:
    message: "Plot interaction profiles from the targets' point of view, inside the exons and 100 nts upstream and downstream."
    input:
        
    output:
    params:
        config["data"]["annotation"]
    script:
        "../scripts/nts_flanking_exons.py"