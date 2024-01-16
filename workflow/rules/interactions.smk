rule sno_RBP_overlap:
    input:
        rules.bedtools_merge_snoGloBe.output
    output:
        "results/interactions/sno_RBP_target_overlap/{sno}_RBP.tsv"
    params:
        gtf = config["path"]["gtf"],
        genome = config["path"]["chrNameLength"],
        preprocessed_ENCODE = "results/interactions/ENCODE"
    conda:
        "../envs/bedtools.yaml"
    message:
        "Calculate p-values for {wildcards.sno}-RBP target overlap interactions."
    shell:
        "bash workflow/scripts/compute_overlaps.sh {input} {params.preprocessed_ENCODE} {params.gtf} {params.genome} {output}"

"""
rule extract_sig_sno_RBP_overlap:
    input:
        in_sno_rbp = expand(rules.sno_rbp_overlap.output,sno=sno_list)
    output:
        temp(os.path.join(config["outpath"],"temp_significant_sno_rbp_target_overlaps.tsv"))
        #out_sno_rbp = temp(os.path.join(config["temp"],"significant_sno_rbp_target_overlaps.tsv")),
        #out_rbp_rbp = temp(os.path.join(config["temp"],"significant_rbp_rbp_target_overlaps.tsv")),
        #out_sno_sno = temp(os.path.join(config["temp"],"significant_sno_sno_target_overlaps.tsv"))
    params:
        config["filters"]["ovlp_p_val_threshold"]
    message:
        "Extract significant snoRNA-RBP target overlap interactions."
    shell:
        "awk \'($3==\"0\") && ($4=={params}) {{print}}\' {input} >> {output}"
        

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


rule filter_STRING:
    input:
        interactions = config['path']['STRING_interactions'],
        info = config['path']['STRING_info']
    output:
        interactions = "results/interactions/STRING/protein_physical_interactions_700.tsv",
        info = "results/interactions/STRING/protein_info.tsv"
    params:
        score_thres = 700
    message:
        "Extract STRING physical binding interactions that are above the combined score threshold."
    shell:
        "echo -e \"protein1\tprotein2\tcombined_score\" > {output.interactions} && "
        "gunzip -c {input.interactions} | awk \'(NR>1 && $3>{params.score_thres}) {{print}}\' >> {output.interactions} && "
        "gunzip -c {input.info} > {output.info}"


rule extract_STRING:
    input:
        interactions = rules.filter_STRING.output.interactions,
        info = rules.filter_STRING.output.info
    output:
        "results/interactions/STRING/physical_RBP_RBP_700.tsv"
    params:
        RBP_list = config['path']['RBP_list']
    conda:
        "../envs/python.yaml"
    message:
        "Extract RBP-RBP physical binding interactions of interest."
    script:
        "../scripts/STRING.py"