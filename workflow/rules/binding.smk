rule get_rbp_transcripts:
    """ Get transcripts for all RBPs in the list """
    input:
        annotation = config["data"]["annotation"],
        rbp = config["nodes"]["rbp_list"],
        sno = config["nodes"]["snoRNA_list"]
    output:
        temp(os.path.join(config["outpath"],"rbp_transcripts.tsv"))
    shell:
        'python scripts/get_transcripts.py rbp {input.annotation} {input.rbp} {input.sno} {output}'

rule get_sno_transcripts:
    message: "Get transcripts for all snoRNAs in the list."
    input:
        annotation = config["data"]["annotation"],
        rbp = config["nodes"]["rbp_list"],
        sno = config["nodes"]["snoRNA_list"]
    output:
        temp(os.path.join(config["outpath"],"sno_transcripts.tsv"))
    shell:
        'python scripts/get_transcripts.py sno {input.annotation} {input.rbp} {input.sno} {output}' 

rule sno_rbp_transcript_bedtools_intersect:
    """ bedtools intersect snoglobe predictions and RBP transcripts """
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

rule rbp_sno_transcript_bedtools_intersect:
    message: "bedtools intersect ENCODE RBP interactions and snoRNA transcripts."
    input:
        left = os.path.join(config["filtered_data"]["rbp_formatted"],"{rbp}_uniq_regions.bed"), ### FIX BACK TO: rules.rbp_final_sort.output,
        right = rules.get_sno_transcripts.output
    output:
        temp(os.path.join(config["outpath"],"{rbp}_sno_transcript_intersect.tsv"))
    params:
        extra = "-wo -s"
    log:
        os.path.join(config["logs"],"rbp_sno_transcript_bedtools_intersect","{rbp}_intersect.log")
    wrapper:
        "v1.3.2/bio/bedtools/intersect"

rule sno_rbp_transcript_out_file:
    """ Format interactions to fit network data format """
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

rule rbp_sno_transcript_out_file:
    """ Format interactions to fit network data format """
    input:
        expand(rules.rbp_sno_transcript_bedtools_intersect.output,rbp=rbp_list)
    output:
        interactions = os.path.join(config["outpath"],"rbp_bind_to_sno_transcript.tsv"),
        counts = temp(os.path.join(config["temp"],"rbp_bind_to_sno_transcript_count.tsv"))
    params:
        temp_dir = config["temp"],
        temp_file = os.path.join(config["temp"],"temp_rbp_bind_to_sno_transcript.tsv")
    shell:
        "echo -e \"RBP\tsnoRNA\tinteraction\" >> {params.temp_file} && "
        "awk '{{print $4\"\t\"$10\"\trbp_sno_transcript\"}}' {params.temp_dir}*_sno_transcript_intersect.tsv >> {params.temp_file} && "
        "awk -F\'\t\' \'{{sub(/\_.+$/,\"\",$1)}}1\' OFS=\'\t\' {params.temp_file} | uniq > {output.interactions} && rm {params.temp_file} && "
        "int_cnt=$(grep -c ^ {output.interactions}) && int_cnt=$(($int_cnt-1)) && echo -e \"{rule}\t$int_cnt\" >> {output.counts}"