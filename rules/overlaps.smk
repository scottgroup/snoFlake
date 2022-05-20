rule sno_rbp_overlap:
    """ Caculate p-value for each snoRNA-RBP overlapping target interaction. """
    input:
        rules.snoglobe_uniq.output
    output:
        os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals","{sno}_rbp_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["data"]["rbp_formatted"],
        calc_script = "scripts/compute_overlaps.py",
        outdir = os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "cd $SLURM_TMPDIR/ && "
        "python3 scripts/compute_overlaps.py {input} {params.rbp_path} {params.gtf} "
        "{params.genome} {output} $SLURM_CPUS_PER_TASK {params.tmpdir} sno_rbp"

rule sno_sno_overlap:
    """ Caculate p-value for each snoRNA-snoRNA overlapping target interaction. """
    input:
        rules.snoglobe_uniq.output
    output:
        os.path.join(config["outpath"],"sno_sno_overlaps_p_vals","{sno}_sno_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        sno_path = config["data"]["snoglobe_formatted"],
        tmp_outfile = "$SLURM_TMPDIR/{sno}_sno_overlaps.tsv",
        calc_script = "scripts/compute_overlaps.py"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "cp {input} {params.gtf} {params.genome} $SLURM_TMPDIR/ && "
        "cp -r {params.sno_path} $SLURM_TMPDIR/ && "    
        "cd $SLURM_TMPDIR/ && "
        ## get basenames of inputs
        "python3 compute_overlaps.py $snofile $tmp_rbp_path $gtf $genome $tmp_outfile $SLURM_CPUS_PER_TASK $SLURM_TMPDIR sno_rbp;"
        "cp {params.tmp_outfile} {output}"    

rule rbp_rbp_overlap:
    """ Caculate p-value for each RBP-RBP overlapping target interaction. """
    input:
        rules.rbp_final_sort.output
    output:
        os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals","{rbp}_rbp_overlaps.tsv")
    params:
        gtf = config["data"]["annotation"],
        genome = config["data"]["chrLength"],
        rbp_path = config["data"]["rbp_formatted"],
        calc_script = "scripts/compute_overlaps.py",
        outdir = os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "cd $SLURM_TMPDIR/ && "
        "python3 scripts/compute_overlaps.py {input} {params.rbp_path} {params.gtf} "
        "{params.genome} {output} $SLURM_CPUS_PER_TASK {params.tmpdir} rbp"

rule extract_sig_overlap:
    """ Extract target overlap interactions that are above limit of detection. """
    input:
        rules.sno_rbp_overlap.output,
        rules.rbp_rbp_overlap.output,
        rules.sno_sno_overlap.output
    output:

    params:
        config["ovlp_p_val_threshold"]
    shell:
        "awk '($3=="0") && ($4=="100000") {print}' $f >> $sig_p_val"