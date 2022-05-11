rule get_rbp_transcripts:
    """ Get transcripts for all RBPs in the list """
    input:
        annotation = config["data"]["annotation"],
        rbp = config["data"]["rbp_list"],
        sno = config["data"]["snoRNA_list"]
    output:
        temp(os.path.join(config["temp"],"rbp_transcripts.tsv"))
    shell:
        'python scripts/get_transcripts.py rbp {input.annotation} {input.rbp} {input.sno} {output}'

rule get_sno_transcripts:
    """ Get transcripts for all snoRNAs in the list """
    input:
        annotation = config["data"]["annotation"],
        rbp = config["data"]["rbp_list"],
        sno = config["data"]["snoRNA_list"]
    output:
        temp(os.path.join(config["temp"],"sno_transcripts.tsv"))
    shell:
        'python scripts/get_transcripts.py sno {input.annotation} {input.rbp} {input.sno} {output}' 

rule sno_rbp_transcript_bedtools_intersect:
    """ bedtools intersect snoglobe predictions and RBP transcripts """
    input:
        left = rules.snoglobe_uniq.output,
        right = rules.get_rbp_transcripts.output
    output:
        temp(os.path.join(config["temp"],"{sno}_rbp_transcript_intersect.tsv"))
    params:
        extra = "-wo -s"
    log:
        os.path.join(config["logs"],"sno_rbp_transcript_bedtools_intersect","{sno}_intersect.log")
    wrapper:
        "v1.3.2/bio/bedtools/intersect"

rule rbp_sno_transcript_bedtools_intersect:
    """ bedtools intersect ENCODE RBP interactions and snoRNA transcripts """
    input:
        left = rules.rbp_final_sort.output,
        right = rules.get_sno_transcripts.output
    output:
        temp(os.path.join(config["temp"],"{rbp}_sno_transcript_intersect.tsv"))
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
        os.path.join(config["outpath"],"sno_bind_to_rbp_transcript.tsv")
    shell:
        "awk '{{print $4\"\t\"$10\"\tsno_rbp_transcript\"}}' {input} | uniq >> {output}"

rule rbp_sno_transcript_out_file:
    """ Format interactions to fit network data format """
    input:
        expand(rules.rbp_sno_transcript_bedtools_intersect.output,rbp=rbp_list)
    output:
        os.path.join(config["outpath"],"rbp_bind_to_sno_transcript.tsv")
    shell:
        "name=$(basename {input} | sed 's/_uniq_regions.bed//g')"
        "awk -v var=name '{{print var\"\t\"$10\"\trbp_sno_transcript\"}}' {input} | uniq >> {output}"