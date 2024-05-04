rule sno_RBP_overlap:
    # Need to use virtualenv instead of conda for this rule as it requires high amount of I/O operations
    input:
        cleanup = rules.preprocessing_cleanup.output,
        sno = "results/interactions/snoGloBe_HTRRI_{sno_thres}/{sno}.bed"
    output:
        "results/interactions/sno_RBP_target_overlap/sno-{sno_thres}-RBP-{rbp_thres}-ovlp-{ovlp_thres}/{sno}_RBP.tsv"
    params:
        gtf = config["path"]["gtf"],
        genome = config["path"]["chrNameLength"],
        preprocessed_ENCODE = "results/interactions/ENCODE_{rbp_thres}",
        virtualenv = config["path"]["virtualenv"],
        env_req = "workflow/envs/overlap_requirements.txt",
        max_iters = "{ovlp_thres}"
    message:
        "Calculate p-values for {wildcards.sno}-RBP target overlap interactions. (THRESHOLDS: {wildcards.sno_thres}, {wildcards.rbp_thres}, {wildcards.ovlp_thres})"
    shell:
        "module load bedtools && "
        "source {params.virtualenv} && "
        "pip install -r {params.env_req} && "
        "bash workflow/scripts/compute_overlaps.sh {input.sno} {params.preprocessed_ENCODE} {params.gtf} {params.genome} {output} {params.max_iters}"


rule extract_sig_sno_RBP_overlap:
    input:
        expand(rules.sno_RBP_overlap.output,sno=sno_list,allow_missing=True)
    output:
        "results/interactions/sno_RBP_target_overlap/sno-{sno_thres}-RBP-{rbp_thres}-ovlp-{ovlp_thres}/all_sig_sno_RBP_target_overlap.tsv"
    params:
        p_val_thres = rules.sno_RBP_overlap.params.max_iters
    message:
        "Extract significant snoRNA-RBP target overlap interactions. (THRESHOLDS: {wildcards.sno_thres}, {wildcards.rbp_thres}, {wildcards.ovlp_thres})"
    shell:
        "echo -e \"source\ttarget\tinteraction\" > {output} && " 
        "awk \'($3==\"0\") && ($4=={params.p_val_thres}) {{print $1\"\t\"$2\"\tsno_RBP_overlap\"}}\' {input} >> {output}"


rule compute_weights_sno_RBP_overlap:
    input:
        expand(rules.extract_sig_sno_RBP_overlap.output,zip,sno_thres=config["thresholds"]["snoGloBe_HTRRI"],rbp_thres=config["thresholds"]["ENCODE"],ovlp_thres=config["thresholds"]["ovlp"])
    output:
        "results/interactions/sno_RBP_target_overlap/weighted_sig_sno_RBP_target_overlap.tsv"
    conda:
        "../envs/bedtools.yaml"
    message:
        "Add weights to significant snoRNA-RBP target overlap interactions."
    script:
        "../scripts/compute_overlaps_weights.py"
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
        score_thres = config["thresholds"]["STRING"]
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
        score_thres = rules.filter_STRING.params.score_thres
    conda:
        "../envs/python.yaml"
    message:
        "Extract RBP-RBP physical binding interactions of interest."
    script:
        "../scripts/STRING.py"


rule RBP_binds_to_sno:
    input:
        expand(rules.bedtools_merge_ENCODE.output,thres=config["thresholds"]["ENCODE"][0],rbp=rbp_list)
    output:
        "results/interactions/RBP_binds_to_sno/all_RBP_binding_to_sno.tsv"
    params:
        preprocessed_ENCODE = "results/interactions/ENCODE_"+str(config["thresholds"]["ENCODE"][0]),
        sno_annotation = config['path']['snoRNA_list'],
        fasta = config['path']['fasta'],
        p_val_thres = config["thresholds"]["ENCODE"][0]
    conda:
        "../envs/bedtools.yaml"
    log:
        "results/logs/RBP_binds_to_sno.log"
    message:
        "Find ENCODE RBPs that bind to snoRNAs."
    shell:
        "python3 workflow/scripts/RBP_binds_to_sno.py {params.preprocessed_ENCODE} "
        "{params.sno_annotation} {params.fasta} {output} {params.p_val_thres} > {log}"
