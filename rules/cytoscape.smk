rule build_network:
    """ Build snoRNA-RBP interaction network using cytoscape. """
    input:
    output:
    conda:
        "../envs/cytoscape.yaml"
    shell:
        "python scripts/cytoscape.py"