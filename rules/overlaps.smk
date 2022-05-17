rule sno_rbp_overlap:
    """ Caculate p-value for each snoRNA-RBP overlapping target interaction. """
    input:
        rules.snoglobe_uniq.output
    output:
        os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals","{sno}.tsv")
    params:
        tmpdir = config["slurm_tmpdir"],
        cpus_per_task = config["slurm_cpus_per_task"],
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["data"]["rbp_formatted"],
        outdir = os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "python3 scripts/compute_overlaps.py {input} {params.rbp_path} {params.gtf} "
        "{params.genome} {output} {params.cpus_per_task} {params.tmpdir} sno_rbp"
"""
rule sno_sno_overlap:
    \"\"\" Caculate p-value for each snoRNA-snoRNA overlapping target interaction. \"\"\"
    input:
        rules.snoglobe_uniq.output
    output:
        os.path.join(config["outpath"],"sno_sno_overlaps_p_vals","{sno}.tsv")
    params:
        tmpdir = config["slurm_tmpdir"],
        cpus_per_task = config["slurm_cpus_per_task"],
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        sno_path = config["data"]["snoglobe_formatted"],
        outdir = os.path.join(config["outpath"],"sno_sno_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "python3 scripts/compute_overlaps.py {input} {params.sno_path} {params.gtf} "
        "{params.genome} {output} {params.cpus_per_task} {params.tmpdir} sno"

rule rbp_rbp_overlap:
    \"\"\" Caculate p-value for each RBP-RBP overlapping target interaction. \"\"\"
    input:
        rules.rbp_final_sort.output
    output:
        os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals","{rbp}.tsv")
    params:
        tmpdir = config["slurm_tmpdir"],
        cpus_per_task = config["slurm_cpus_per_task"],
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["data"]["rbp_formatted"],
        outdir = os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "python3 scripts/compute_overlaps.py {input} {params.rbp_path} {params.gtf} "
        "{params.genome} {output} {params.cpus_per_task} {params.tmpdir} rbp"
"""