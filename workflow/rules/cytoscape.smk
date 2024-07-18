rule network_data:
    input:
        RBP_binds_to_sno = rules.RBP_binds_to_sno.output,
        STRING = rules.extract_STRING.output,
        sno_RBP_overlap = rules.bh_sno_RBP_overlap.output.bh,
        snoRNA = config["path"]["snoRNA_list"],
        RBP = config["path"]["RBP_list"],
    output:
        nodes = "results/networks/nodes.tsv",
        edges = "results/networks/edges.tsv"
    params:
        type = "data"
    conda:
        "../envs/python.yaml"
    message:
        "Gather node and edge data to construct the snoRNA-RBP interaction network via Cytoscape."
    shell:
        "pip install py4cytoscape && "
        "python3 workflow/scripts/cytoscape.py {params.type} {output.nodes} {output.edges} {input.RBP_binds_to_sno} {input.STRING} "
        "{input.sno_RBP_overlap} {input.snoRNA} {input.RBP}"


rule build_network:
    # Need Cytoscape app to be RUNNING
    input:
        nodes = rules.network_data.output.nodes,
        edges = rules.network_data.output.edges
    output:
        network = "results/networks/snoFlake.cys"
    params:
        type = "build"
    conda:
        "../envs/python.yaml"
    message:
        "Build the snoRNA-RBP interaction network via Cytoscape."
    shell:
        "python3 workflow/scripts/cytoscape.py {params.type} {input.nodes} {input.edges} {output.network}"