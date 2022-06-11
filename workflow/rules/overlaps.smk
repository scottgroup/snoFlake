
rule sno_rbp_overlap:
    message: "Calculate p-value for each snoRNA-RBP overlapping target interaction."
    input:
        rules.sno_final_sort.output
    output:
        os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals","{sno}_rbp_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["filtered_data"]["rbp_formatted"],
        sno = "{sno}",
        outdir = os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "echo {params.sno} && "
        "bash workflow/scripts/compute_overlaps.sh {input} {params.rbp_path} {params.gtf} {params.genome} {output} sno_rbp && "
        "echo \'Done\'"
"""
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
        "bash scripts/compute_overlaps.sh {input} {params.sno_path} {params.gtf} {params.genome} {output} sno && "
        "echo \'Done\'"  

rule rbp_rbp_overlap:
    message: "Caculate p-value for each RBP-RBP overlapping target interaction."
    input:
        os.path.join(config["data"]["rbp_formatted"],"NSUN2_uniq_regions.bed")
        #rules.rbp_final_sort.output
    output:
        os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals","NSUN2_rbp_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["data"]["rbp_formatted"],
        rbp = "NSUN2",
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
        in_sno_rbp = expand(rules.sno_rbp_overlap.output,sno=sno_list)
        #in_rbp_rbp = expand(rules.rbp_rbp_overlap.output,rbp=rbp_list),
        #in_sno_sno = expand(rules.sno_sno_overlap.output,sno=sno_list)
    output:
        temp(os.path.join(config["outpath"],"temp_significant_sno_rbp_target_overlaps.tsv"))
        #out_sno_rbp = temp(os.path.join(config["temp"],"significant_sno_rbp_target_overlaps.tsv")),
        #out_rbp_rbp = temp(os.path.join(config["temp"],"significant_rbp_rbp_target_overlaps.tsv")),
        #out_sno_sno = temp(os.path.join(config["temp"],"significant_sno_sno_target_overlaps.tsv"))
    params:
        config["filters"]["ovlp_p_val_threshold"]
    shell:
        "awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input} >> {output}"
        #"awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input.in_sno_rbp} >> {output.out_sno_rbp}; "
        #"awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input.in_rbp_rbp} >> {output.out_rbp_rbp}; "
        #"awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input.in_sno_sno} >> {output.out_sno_sno}; "

rule format_outfile:
    message: "Format output files to meet network configurations."
    input:
        rules.extract_sig_overlaps.output
        #in_sno_rbp = rules.extract_sig_overlaps.output.out_sno_rbp,
        #in_rbp_rbp = rules.extract_sig_overlaps.output.out_rbp_rbp,
        #in_sno_sno = rules.extract_sig_overlaps.output.out_sno_sno
    output:
        os.path.join(config["outpath"],"significant_sno_rbp_target_overlaps.tsv")
        #out_sno_rbp = os.path.join(config["outpath"],"significant_sno_rbp_target_overlaps.tsv"),
        #out_rbp_rbp = os.path.join(config["outpath"],"significant_rbp_rbp_target_overlaps.tsv"),
        #out_sno_sno = os.path.join(config["outpath"],"significant_sno_sno_target_overlaps.tsv")
    shell:
        "echo -e \"source\ttarget\tinteraction\" >> {output} && " 
        "sed 's/.bed//g' {input} | awk \'{{print $1\"\t\"$2\"\t\"\"sno_rbp_target_overlap\"}}\' >> {output}; "
        #"echo -e \"snoRNA\tRBP\tinteraction\" >> {output.out_sno_rbp} && "
        #"awk \'{{print $1\"\t\"$2\"\t\"\"sno_rbp_target_overlap\"}}\' {input.in_sno_rbp} >> {output.out_sno_rbp}; "
        #"echo -e \"RBP1\tRBP2\tinteraction\" >> {output.out_rbp_rbp} && "
        #"awk \'{{print $1\"\t\"$2\"\t\"\"rbp_rbp_target_overlap\"}}\' {input.in_rbp_rbp} >> {output.out_rbp_rbp}; "
        #"echo -e \"snoRNA1\tsnoRNA2\tinteraction\" >> {output.out_sno_sno} && "
        #"awk \'{{print $1\"\t\"$2\"\t\"\"sno_sno_target_overlap\"}}\' {input.in_sno_sno} >> {output.out_sno_sno}; "
"""
rule get_targets:
    message: "Get overalpping targets."
    input:
        rule.format_outfile.output
    output:
        os.path.join(config["outpath"],"sno_rbp_overlaps_targets")
    
    shell:
        "awk '(NR>1) {print $1\"\t\"$2}' {input} | while read line; do"
        "done"

def sig_pairs(file):
    source = pd.read_table(file).source.values.tolist()
    target = pd.read_table(file).target.values.tolist()
    return expand("histone_{histone}/bin_size_{bin_size}/result.csv", zip,
       source_node=source, target_node=target)

rule get_overlapping_targets_bedtools_intersect:
    message: "bedtools intersect to get overlapping targets."
    input:
        left = 
        right = 
    output:
    log:
        os.path.join(config["logs"],"get_overlapping_targets_bedtools_intersect","{rbp}.log")
"""