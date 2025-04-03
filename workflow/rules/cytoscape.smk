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


rule network_motifs:
    # Need Cytoscape app to be RUNNING
    input:
        edges = rules.network_data.output.edges
    output:
        network_motifs = "results/networks/network_motifs.tsv"
    conda:
        "../envs/python.yaml"
    message:
        "Search for double-edged network motifs."
    script:
        "../scripts/network_motifs.py"


rule build_network:
    # Need Cytoscape app to be RUNNING
    input:
        nodes = rules.network_data.output.nodes,
        edges = rules.network_data.output.edges,
        motifs = rules.network_motifs.output.network_motifs
    output:
        network = "results/networks/snoFlake.cys"
    params:
        type = "build"
    conda:
        "../envs/python.yaml"
    message:
        "Build the snoRNA-RBP interaction network via Cytoscape."
    shell:
        "python3 workflow/scripts/cytoscape.py {params.type} {input.nodes} {input.edges} {input.motifs} {output.network}"


rule extra_node_annotations:
    input:
        nodes = rules.network_data.output.nodes,
        encore_matrix = config['path']['ENCORE'],
        rRNA_targets = config['path']['rRNA_targets'],
        snRNA_targets = config['path']['snRNA_targets'],
        snodb = config['path']['snoDB']
    output:
        RBP_anno = "results/networks/RBP_function_localization.tsv",
        sno_anno = "results/networks/snoRNA_targets.tsv"
    conda:
        "../envs/python.yaml"
    message:
        "Add RBP function and localization data as well as snoRNA target information."
    script:
        "../scripts/node_annotation.py"


rule sno_RBP_overlap_targets_annotations:
    input:
        targets = expand(rules.sno_RBP_overlap_targets.output,sno=sno_list),
        nodes = rules.network_data.output.nodes,
        edges = rules.network_data.output.edges,
        exp_pc_genes = config['path']['pc_list']
    output:
        nodes = "results/networks/nodes_targets.tsv",
        edges = "results/networks/edges_targets.tsv"
    params:
        target_path = "results/interactions/sno_RBP_target_overlap/targets"
    conda:
        "../envs/python.yaml"
    message:
        "Generate node and edge file for significant snoRNA-RBP overlapping targets."
    script:
        "../scripts/cytoscape_targets.py"