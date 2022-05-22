"""
rule sno_rbp_overlap:
    message: "Caculate p-value for each snoRNA-RBP overlapping target interaction."
    input:
        rules.snoglobe_uniq.output
    output:
        os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals","{sno}_rbp_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["data"]["rbp_formatted"],
        sno = "{sno}",
        outdir = os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "echo {params.sno} && "
        "mkdir -p {params.outdir} && "
        "scripts/compute_overlaps.sh {input} {params.rbp_path} {params.gtf} {params.genome} {output} sno_rbp && "
        "echo \'Done\'"

rule sno_sno_overlap:
    message: "Caculate p-value for each snoRNA-snoRNA overlapping target interaction."
    input:
        rules.snoglobe_uniq.output
    output:
        os.path.join(config["outpath"],"sno_sno_overlaps_p_vals","{sno}_sno_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        sno_path = config["data"]["snoglobe_formatted"],
        sno = "{sno}",
        outdir = os.path.join(config["outpath"],"sno_sno_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "echo {params.sno} && "
        "mkdir -p {params.outdir} && "
        "scripts/compute_overlaps.sh {input} {params.sno_path} {params.gtf} {params.genome} {output} sno && "
        "echo \'Done\'"  

rule rbp_rbp_overlap:
    message: "Caculate p-value for each RBP-RBP overlapping target interaction."
    input:
        rules.rbp_final_sort.output
    output:
        os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals","{rbp}_rbp_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["data"]["rbp_formatted"],
        rbp = "{rbp}",
        outdir = os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "echo {params.rbp} && "
        "mkdir -p {params.outdir} && "
        "scripts/compute_overlaps.sh {input} {params.rbp_path} {params.gtf} {params.genome} {output} rbp && "
        "echo \'Done\'"  
"""
rule extract_sig_overlaps:
    message: "Extract target overlap interactions that are above limit of detection."
    input:
        expand(os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals","{sno}_rbp_overlap.tsv"),sno=sno_list)
        #in_sno_rbp = expand(rules.sno_rbp_overlap.output,sno=sno_list),
        #in_rbp_rbp = expand(rules.rbp_rbp_overlap.output,rbp=rbp_list),
        #in_sno_sno = expand(rules.sno_sno_overlap.output,sno=sno_list)
    output:
        temp(os.path.join(config["temp"],"significant_snoglobe_sno_rbp_target_overlaps.tsv"))
        #out_sno_rbp = temp(os.path.join(config["temp"],"significant_sno_rbp_target_overlaps.tsv")),
        #out_rbp_rbp = temp(os.path.join(config["temp"],"significant_rbp_rbp_target_overlaps.tsv")),
        #out_sno_sno = temp(os.path.join(config["temp"],"significant_sno_sno_target_overlaps.tsv"))
    params:
        config["ovlp_p_val_threshold"]
    shell:
        "awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input} >> {output}"
        #"awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input.in_sno_rbp} >> {output.out_sno_rbp}; "
        #"awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input.in_rbp_rbp} >> {output.out_rbp_rbp}; "
        #"awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input.in_sno_sno} >> {output.out_sno_sno}; "

rule format_outfile:
    message: "Format output files to meet network requirements."
    input:
        rules.extract_sig_overlaps.output
        #in_sno_rbp = rules.extract_sig_overlaps.output.out_sno_rbp,
        #in_rbp_rbp = rules.extract_sig_overlaps.output.out_rbp_rbp,
        #in_sno_sno = rules.extract_sig_overlaps.output.out_sno_sno
    output:
        os.path.join(config["outpath"],"significant_snoglobe_sno_rbp_target_overlaps.tsv")
        #out_sno_rbp = os.path.join(config["outpath"],"significant_sno_rbp_target_overlaps.tsv"),
        #out_rbp_rbp = os.path.join(config["outpath"],"significant_rbp_rbp_target_overlaps.tsv"),
        #out_sno_sno = os.path.join(config["outpath"],"significant_sno_sno_target_overlaps.tsv")
    shell:
        "echo -e \"snoRNA\tRBP\tinteraction\" >> {output} && "
        "awk \'{{print $1\"\t\"$2\"\t\"\"snoglobe_sno_rbp_target_overlap\"}}\' {input} >> {output}; "
        #"echo -e \"snoRNA\tRBP\tinteraction\" >> {output.out_sno_rbp} && "
        #"awk \'{{print $1\"\t\"$2\"\t\"\"sno_rbp_target_overlap\"}}\' {input.in_sno_rbp} >> {output.out_sno_rbp}; "
        #"echo -e \"RBP1\tRBP2\tinteraction\" >> {output.out_rbp_rbp} && "
        #"awk \'{{print $1\"\t\"$2\"\t\"\"rbp_rbp_target_overlap\"}}\' {input.in_rbp_rbp} >> {output.out_rbp_rbp}; "
        #"echo -e \"snoRNA1\tsnoRNA2\tinteraction\" >> {output.out_sno_sno} && "
        #"awk \'{{print $1\"\t\"$2\"\t\"\"sno_sno_target_overlap\"}}\' {input.in_sno_sno} >> {output.out_sno_sno}; "