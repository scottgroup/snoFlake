rule build_network:
    message: "Build snoRNA-RBP interaction network using cytoscape."
    input:
    output:
    conda:
        "../envs/cytoscape.yaml"
    script:
        "../scripts/cytoscape.py"