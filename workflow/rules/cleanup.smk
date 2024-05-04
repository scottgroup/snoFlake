rule preprocessing_cleanup:
    input:
        ENCODE = expand(rules.bedtools_merge_ENCODE.output,thres=config["thresholds"]["ENCODE"],rbp=rbp_list),
        snoGloBe = expand(rules.bedtools_merge_snoGloBe_HTRRI.output,thres=config["thresholds"]["snoGloBe_HTRRI"],sno=sno_list)
    output:
        "results/preprocessing_status.txt"
    message:
        "Remove temp files created during preprocessing."
    shell:
        "find -type f -path '*/ENCODE_*/*_filtered_merged.bed' -delete && "
	    "echo -e \"Preprocessing of snoRNA and RBP interactions complete.\" > {output} && "
        "echo -e \"Computing network interactions.\" >> {output}"

"""
rule interactions_cleanup:
    input:
        STRING = rules.extract_STRING.output
    output:
        "results/network_status.txt"
    message:
        "Remove temp files created during network interaction computation."
    shell:
"""