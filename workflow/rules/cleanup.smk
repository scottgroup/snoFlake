rule cleanup_files:
    input:
        ENCODE = expand(rules.bedtools_merge_ENCODE.output,rbp=rbp_list),
        snoGloBe = expand(rules.bedtools_merge_snoGloBe_HTRRI.output,sno=sno_list)
    output:
        "results/preprocessing_status.txt"
    message:
        "Remove temp files."
    shell:
        "if test -f results/interactions/ENCODE/*_filtered_merged.bed; then rm results/interactions/ENCODE/*_filtered_merged.bed; fi  && "
	"echo -e \"Preprocessing of snoRNA and RBP interactions complete.\" > {output} && "
        "echo -e \"Computing network interactions.\" >> {output}"
