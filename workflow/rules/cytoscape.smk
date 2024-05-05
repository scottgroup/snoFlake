rule build_network:
    input:
        RBP_binds_to_sno = rules.RBP_binds_to_sno.output,
        STRING = rules.extract_STRING.output,
        sno_RBP_overlap = rules.compute_weights_sno_RBP_overlap.output,
        snoRNA = config["path"]["snoRNA_list"],
        RBP = config["path"]["RBP_list"],
        pc_gene = config["path"]["exp_pc_genes"]
    output:
        network = "results/networks/sno_RBP_network.cys",
        nodes = "results/networks/nodes.tsv",
        edges = "results/networks/edges.tsv"
    params:
        virtualenv = config["path"]["virtualenv"],
        env_req = "workflow/envs/overlap_requirements.txt"
    message:
        "Build the snoRNA-RBP interaction network via Cytoscape."
    shell:
        "source {params.virtualenv} && "
        "pip install -r {params.env_req} && "
        "python3 workflow/scripts/cytoscape.py {input.snoRNA} {input.RBP} {input.pc_gene} {input.STRING} {input.sno_RBP_overlap} "
        "{input.RBP_binds_to_sno} {output.nodes} {output.edges} {output.network}"