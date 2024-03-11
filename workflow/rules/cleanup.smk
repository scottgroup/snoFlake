rule cleanup_files:
    input:
        ENCODE = expand(rules.bedtools_merge_ENCODE.output,rbp=rbp_list),
        snoGloBe = expand(rules.bedtools_merge_snoGloBe_HTRRI.output,sno=sno_list)
    output:
        "results/preprocessing_status.txt"
    message:
        "Remove temp files."
    shell:
        "rm results/interactions/ENCODE/*_filtered_merged.bed && "
        "rm results/interactions/snoGloBe_HTRRI/*_tmp.bed && "
        "rm results/interactions/snoGloBe_HTRRI/*.sorted.bed && "
	    "echo -e \"Preprocessing of snoRNA and RBP interactions complete.\" > {output} && "
        "echo -e \"Computing network interactions.\" >> {output}"
