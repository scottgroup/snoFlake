rule build_network:
    input:
        RBP_binds_to_sno = rules.RBP_binds_to_sno.output,
        STRING = rules.extract_STRING.output,
        sno_RBP_overlap = rules.bh_sno_RBP_overlap.output.bh,
        snoRNA = config["path"]["snoRNA_list"],
        RBP = config["path"]["RBP_list"],
    output:
        nodes = "results/networks/nodes.tsv",
        edges = "results/networks/edges.tsv"
    conda:
        "../envs/python.yaml"
    message:
        "Build the snoRNA-RBP interaction network via Cytoscape."
    script:
        "../scripts/cytoscape.py"