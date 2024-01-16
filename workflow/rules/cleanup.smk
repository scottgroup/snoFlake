rule cleanup_files:
    output:
        "results/summary.txt"
    shell:
        "rm results/interactions/ENCODE/*_filtered_merged.bed && "
        "rm results/interactions/snoGloBe/*_tmp.bed && "
        "touch {output}"