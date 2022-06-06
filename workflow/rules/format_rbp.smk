rule filter_rbp:
    message: "Filter raw ENCODE RBP data by p-value and by fold-enrichment and merge replicates."
    input:
        expand(os.path.join(config["data"]["rbp_ENCODE"],"{{encode}}","replicate{num}.bed.gz"),num=["1","2"])
    output:
        temp(os.path.join(config["filtered_data"]["rbp_formatted"],"{encode}_filtered.bed"))
    params:
        p_thres = config["filters"]["ENCODE_p_threshold"], # specified p-value threshold in config.yaml
        sig_val_thres = config["filters"]["ENCODE_sig_val_threshold"] # specified signal value threshold in config.yaml
    shell:
        "gunzip -c {input} | sort -k8 -n | awk -F'\t' '$8>={params.p_thres} && $7>={params.sig_val_thres}' >> {output}"

rule create_empty_merge_file:
    message: "Create empty merge files to automatically delete them after run."
    output:
        temp(os.path.join(config["filtered_data"]["rbp_formatted"],"{rbp}_empty.bed"))
    shell:
        "touch {output}"

rule merge_cell_rbp:
    message: "Merge by cell line for each RBP."
    input:
        expand(rules.filter_rbp.output, encode=ENCODE),
        empty = rules.create_empty_merge_file.output
    output:
        temp(os.path.join(config["filtered_data"]["rbp_formatted"],"{rbp}_all_merge.bed"))
    params:
        config["filtered_data"]["rbp_formatted"]
    shell:
        "ls {params}*_filtered.bed | grep '{wildcards.rbp}_' | while read line; "
        "do "
            "cat ${{line}} >> {output}; "
        "done"

rule sort_rbp:
    message: "Sort RBP data for bedtools merge."
    input:
        rules.merge_cell_rbp.output
    output:
        temp(os.path.join(config["filtered_data"]["rbp_formatted"],"{rbp}_sort.bed"))
    shell:
        "cut -f 1-6 {input} | sort -k1,1 -k2,2n > {output}"

rule bedtools_merge_rbp: 
    message: "bedtools merge sorted RBP data."
    input:
        rules.sort_rbp.output
    output:
        temp(os.path.join(config["filtered_data"]["rbp_formatted"],"{rbp}_bedtools_merge.bed"))
    params:
        extra="-s -c 4,5,6 -o distinct,min,distinct"
    log:
        os.path.join(config["logs"],"bedtools_merge_rbp","{rbp}.log")
    wrapper:
        "v1.3.2/bio/bedtools/merge"

rule rbp_final_sort:
    message: "Remove redundancy in bed file if present."
    input:
        rules.bedtools_merge_rbp.output
    output:
        os.path.join(config["filtered_data"]["rbp_formatted"],"{rbp,[A-Za-z0-9]+}.bed")
    shell:
        "sort -k1,1 -k2,3n -u {input} > {output}"