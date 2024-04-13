"""
rule sno_RBP_overlap:
    # Need to use virtualenv instead of conda for this rule as it requires high amount of I/O operations
    input:
        sno = rules.bedtools_merge_snoGloBe_HTRRI.output,
        cleanup = rules.cleanup_files.output
    output:
        "results/interactions/sno_RBP_target_overlap/{sno}_RBP.tsv"
    params:
        gtf = config["path"]["gtf"],
        genome = config["path"]["chrNameLength"],
        preprocessed_ENCODE = "results/interactions/ENCODE",
        virtualenv = config["path"]["virtualenv"],
        env_req = "workflow/envs/overlap_requirements.txt"
    message:
        "Calculate p-values for {wildcards.sno}-RBP target overlap interactions."
    shell:
        "module load bedtools && "
        "source {params.virtualenv} && "
        "pip install -r {params.env_req} && "
        "bash workflow/scripts/compute_overlaps.sh {input.sno} {params.preprocessed_ENCODE} {params.gtf} {params.genome} {output}"
"""

rule extract_sig_sno_RBP_overlap:
    input:
        #expand(rules.sno_RBP_overlap.output,sno=sno_list)
        expand("results/interactions/sno_RBP_target_overlap/{sno}_RBP.tsv",sno=sno_list)
    output:
        "results/interactions/sno_RBP_target_overlap/all_sig_sno_RBP_target_overlap.tsv"
    params:
        p_val_thres = 100000
    message:
        "Extract significant snoRNA-RBP target overlap interactions."
    shell:
        "echo -e \"source\ttarget\tinteraction\" > {output} && " 
        "awk \'($3==\"0\") && ($4=={params.p_val_thres}) {{print $1\"\t\"$2\"\tsno_RBP_overlap\"}}\' {input} >> {output}"

"""
rule exp_overlapping_targets:
    input:
        sno_interaction = "results/interactions/snoGloBe/ENSG00000277194.bed",
        sno_RBP_overlap = "results/interactions/sno_RBP_target_overlap/ENSG00000277194_RBP.tsv"
    output:
        "results/overlapping_targets/{sno}/summary.txt"
    params:
        p_val_thres = 100000,
        exp_pc_genes = config['path']['exp_pc_genes'],
        preprocessed_ENCODE = "results/interactions/ENCODE",
        outdir = "results/overlapping_targets/{sno}"
    conda:
        "../envs/bedtools.yaml"
    message:
        "Find overlapping targets of {wildcards.sno}-RBP that are expressed."
    script:
        "../scripts/exp_overlapping_targets.py"
"""

rule filter_STRING:
    input:
        interactions = config['path']['STRING_interactions'],
        info = config['path']['STRING_info']
    output:
        interactions = "results/interactions/STRING/protein_physical_interactions.tsv",
        info = "results/interactions/STRING/protein_info.tsv"
    params:
        score_thres = 700
    message:
        "Extract STRING physical binding interactions that are above the combined score threshold."
    shell:
        "echo -e \"protein1\tprotein2\tcombined_score\" > {output.interactions} && "
        "gunzip -c {input.interactions} | awk \'(NR>1 && $3>={params.score_thres}) {{print}}\' >> {output.interactions} && "
        "gunzip -c {input.info} > {output.info}"


rule extract_STRING:
    input:
        interactions = rules.filter_STRING.output.interactions,
        info = rules.filter_STRING.output.info
    output:
        "results/interactions/STRING/physical_RBP_RBP.tsv"
    params:
        RBP_list = config['path']['RBP_list'],
        score_thres = 700
    conda:
        "../envs/python.yaml"
    message:
        "Extract RBP-RBP physical binding interactions of interest."
    script:
        "../scripts/STRING.py"


rule RBP_binds_to_sno:
    input:
        expand(rules.bedtools_merge_ENCODE.output,rbp=rbp_list)
    output:
        "results/interactions/RBP_binds_to_sno/all_RBP_binding_to_sno.tsv"
    params:
        preprocessed_ENCODE = "results/interactions/ENCODE",
        sno_annotation = config['path']['snoRNA_list']
    conda:
        "../envs/bedtools.yaml"
    message:
        "Find ENCODE RBPs that bind to snoRNAs."
    script:
        "../scripts/RBP_binds_to_sno.py"
