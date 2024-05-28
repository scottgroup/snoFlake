rule sno_RBP_overlap:
    input:
        sno = rules.bedtools_merge_snoGloBe_HTRRI.output,
        rbp = expand(rules.bedtools_merge_ENCODE.output, rbp=rbp_list)
    output:
        "results/interactions/sno_RBP_target_overlap/{sno}_RBP.tsv"
    params:
        ENCODE_dir = "results/interactions/ENCODE",
        genome = config["path"]["genome"],
        merged = "results/interactions/sno_RBP_target_overlap/merged_sno_RBP_overlap.tsv"
    conda:
        "../envs/bedtools.yaml"
    message:
        "Calculate p-values for {wildcards.sno}-RBP target overlap interactions."
    shell:
        "echo -e \"left\tright\ttwo-tail\tratio\tsnoRNA\tRBP\" > {output} && "
        "ls {params.ENCODE_dir} | while read line; do "
        "bedtools fisher -a {input.sno} -b {params.ENCODE_dir}/$line -g {params.genome} -s | awk \'NF && $1!~/^#/\' | "
        "awk -v rbp=\"${{line%.*}}\" -v sno={wildcards.sno} \'(NR>1) {{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"sno\"\t\"rbp}}\' >> {output}; "
        "done; "
        "tail -n +2 {output} >> {params.merged}"


rule bh_sno_RBP_overlap:
    input:
        all = expand(rules.sno_RBP_overlap.output, sno=sno_list)
    output:
        bh = "results/interactions/sno_RBP_target_overlap/bh_sno_RBP_overlap.tsv"
    params:
        merged = rules.sno_RBP_overlap.params.merged,
        fdr = config['thresholds']['sno_RBP_overlap_fdr']
    conda:
        "../envs/bedtools.yaml"
    message:
        "Apply Benjamini–Hochberg procedure to find the significance threshold for sno_RBP_overlap p-values."
    script:
        "../scripts/bh_sno_RBP_overlap.py"


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
        "Extract STRING physical binding interactions that are above the combined score threshold {params.score_thres}."
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
        RBP_list = config['path']['RBP_list']
    conda:
        "../envs/bedtools.yaml"
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
        sno_annotation = config['path']['snoRNA_list'],
        fasta = config['path']['fasta']
    conda:
        "../envs/bedtools.yaml"
    log:
        "results/logs/RBP_binds_to_sno.log"
    message:
        "Find ENCODE RBPs that bind to snoRNAs."
    shell:
        "python3 workflow/scripts/RBP_binds_to_sno.py {params.preprocessed_ENCODE} "
        "{params.sno_annotation} {params.fasta} {output} > {log}"
