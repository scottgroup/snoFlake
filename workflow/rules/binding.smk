"""
rule get_rbp_transcripts:
    message: "Get transcripts for all RBPs in the list."
    input:
        annotation = config["data"]["annotation"],
        rbp = config["nodes"]["rbp_list"],
        sno = config["nodes"]["snoRNA_list"]
    output:
        temp(os.path.join(config["outpath"],"rbp_transcripts.tsv"))
    params:
        "rbp"
    script:
        "../scripts/get_transcripts.py"
"""
#### sno --> RBP transcript: extract directly from merge_snoglobe_htrri output and sort | uniq

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

"""
rule sno_rbp_transcript_bedtools_intersect:
    message: "bedtools intersect snoGloBe predictions and RBP transcripts."
    input:
        left = os.path.join(config["filtered_data"]["snoglobe_formatted"],"{sno}_uniq_regions.bed"), ### FIX BACK TO: rules.snoglobe_uniq.output,
        right = rules.get_rbp_transcripts.output
    output:
        temp(os.path.join(config["outpath"],"{sno}_rbp_transcript_intersect.tsv"))
    params:
        extra = "-wo -s"
    log:
        os.path.join(config["logs"],"sno_rbp_transcript_bedtools_intersect","{sno}_intersect.log")
    wrapper:
        "v1.3.2/bio/bedtools/intersect"
"""
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
"""
rule sno_rbp_transcript_out_file:
    message: "Format interactions to fit network data format."
    input:
        expand(rules.sno_rbp_transcript_bedtools_intersect.output,sno=sno_list)
    output:
        interactions = os.path.join(config["outpath"],"snoglobe_sno_bind_to_rbp_transcript.tsv"),
        counts = temp(os.path.join(config["temp"],"snoglobe_sno_bind_to_rbp_transcript_count.tsv"))
    params:
        config["temp"]
    shell:
        "echo -e \"snoRNA\tRBP\tinteraction\" >> {output.interactions} && "
        "awk '{{print $4\"\t\"$10\"\tsno_rbp_transcript\"}}' {params}*_rbp_transcript_intersect.tsv | uniq >> {output.interactions} && "
        "int_cnt=$(grep -c ^ {output.interactions}) && int_cnt=$(($int_cnt-1)) && echo -e \"{rule}\t$int_cnt\" >> {output.counts}"
"""
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