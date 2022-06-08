rule build_network:
    message: "Build snoRNA-RBP interaction network using cytoscape."
    input:
        sno = rules.rank_sno.output,
        rbp = rules.rank_rbp.output,
        rbp_bind_to_sno = rules.rbp_sno_transcript_out_file.output
    output:
    conda:
        "../envs/cytoscape.yaml"
    script:
        "../scripts/cytoscape.py"