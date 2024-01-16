rule cleanup_files:
    input:
        ENCODE = expand(rules.bedtools_merge_ENCODE.output,rbp=rbp_list),
        snoGloBe = expand(rules.bedtools_merge_snoGloBe.output,sno=sno_list)
    output:
        "results/summary.txt"
    message:
        "Remove temp files."
    shell:
        "rm results/interactions/ENCODE/*_filtered_merged.bed && "
        "rm results/interactions/snoGloBe/*_tmp.bed && "
        "echo -e \"Preprocessing of snoRNA and RBP interactions complete.\" > {output}"