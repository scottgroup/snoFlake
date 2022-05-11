rule sno_rbp_overlap:
    input:
        snoglobe = rules.snoglobe_uniq.output,
        gtf = config["data"]["annotation"],
        chr_length = 
    output:
        os.path.join(config["outpath"],"{sno}_uniq_regions.bed")
    params:
        rbp_path = config["data"]["rbp_formatted"]
    shell:
        "module load bedtools/2.30.0 && "
        "module load python/3.8.10 && "

rule sno_sno_overlap:

rule rbp_rbp_overlap: