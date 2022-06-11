rule get_sno_transcripts:
    message: "Get transcripts for all snoRNAs in the list."
    input:
        annotation = config["data"]["annotation"],
        rbp = config["nodes"]["rbp_list"],
        sno = config["nodes"]["snoRNA_list"]
    output:
        temp(os.path.join(config["outpath"],"sno_transcripts.tsv"))
    params:
        "sno"
    script:
        "../scripts/get_transcripts.py"

rule rbp_sno_transcript_bedtools_intersect:
    message: "bedtools intersect ENCODE RBP interactions and snoRNA transcripts."
    input:
        left = rules.rbp_final_sort.output,
        right = rules.get_sno_transcripts.output
    output:
        temp(os.path.join(config["outpath"],"{rbp}_sno_transcript_intersect.tsv"))
    params:
        extra = "-wo -s"
    log:
        os.path.join(config["logs"],"rbp_sno_transcript_bedtools_intersect","{rbp}.log")
    wrapper:
        "v1.3.2/bio/bedtools/intersect"

rule rbp_sno_transcript_out_file:
    message: "Format rbp --> sno transcript interactions to fit network data format."
    input:
        expand(rules.rbp_sno_transcript_bedtools_intersect.output,rbp=rbp_list)
    output:
        os.path.join(config["outpath"],"rbp_bind_to_sno_transcript.tsv")
    params:
        outdir = config["outpath"],
        temp_file = os.path.join(config["outpath"],"temp_rbp_bind_to_sno_transcript.tsv")
    shell:
        "echo -e \"source\ttarget\tinteraction\" >> {params.temp_file} && "
        "awk '{{print $4\"\t\"$10\"\trbp_sno_transcript\"}}' {params.outdir}*_sno_transcript_intersect.tsv >> {params.temp_file} && "
        "awk -F\'\t\' \'{{sub(/\_.+$/,\"\",$1)}}1\' OFS=\'\t\' {params.temp_file} | uniq > {output} && rm {params.temp_file}"

rule sno_pc_targets:
    message: "Get all protein coding targets of snoRNAs."
    input:
        expand(rules.merge_snoglobe_htrri.output,sno=sno_list)
    output:
        temp(os.path.join(config["outpath"],"sno_pc_targets.tsv"))
    shell:
        "awk '{{print $4\"\t\"$12\"\tsno_pc_targets\"}}' {input} | sort -u >> {output}"

rule filter_targets:
    message: "Only keep targets that belong to the RBP list."
    input:
        rules.sno_pc_targets.output
    output:
        os.path.join(config["outpath"],"sno_bind_to_rbp_transcript.tsv")
    params:
        config["nodes"]["rbp_list"]
    shell:
        "echo -e \"source\ttarget\tinteraction\" > {output} && "
        "sed -i \'s/sno_pc_targets/sno_bind_to_rbp_transcript/g\' {input} && "
        "awk \'{{print $1\",\"$2}}\' {params} | while IFS=',' read -r id name; do "
        "sed -i \"s/$id/$name/g\" {input}; done && "
        "awk -F\'\t\' \'$2 !~ \"^(ENSG00000|trna|NR_)\"\' {input} | sed \"s/\,ENSG00000[^\t]*//g\" | sort -u >> {output}"