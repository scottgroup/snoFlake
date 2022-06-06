rule filter_rbp:
    message: "Filter raw ENCODE RBP data by p-value and by fold-enrichment and merge replicates."
    input:
        rep1 = os.path.join(config["data"]["rbp_ENCODE"],"{encode}","replicate1.bed.gz"),
        rep2 = os.path.join(config["data"]["rbp_ENCODE"],"{encode}","replicate2.bed.gz")
    output:
        temp(os.path.join(config["temp"],"{encode}_filtered.bed"))
    params:
        temp_dir = config["temp"],
        p_thres = config["ENCODE_p_threshold"], # specified p-value threshold in config.json
        sig_val_thres = config["ENCODE_sig_val_threshold"] # specified signal value threshold in config.json
    shell:
        "mkdir -p {params.temp_dir} && "
        "gunzip -c {input.rep1} | sort -k8 -n | awk -F'\t' '$8>={params.p_thres} && $7>={params.sig_val_thres}' >> {output} && "
        "gunzip -c {input.rep2} | sort -k8 -n | awk -F'\t' '$8>={params.p_thres} && $7>={params.sig_val_thres}' >> {output}"

rule create_empty_merge_file: # empty temp merge files to automatically delete them after run
    output:
        temp(os.path.join(config["temp"],"{rbp}_empty.bed"))
    shell:
        "touch {output}"

rule merge_cell_rbp: # merge by cell line for each RBP
    input:
        expand(rules.filter_rbp.output, encode=ENCODE),
        empty = rules.create_empty_merge_file.output
    output:
        temp(os.path.join(config["temp"],"{rbp}_all_merge.bed"))
    params:
        config["temp"]
    shell:
        "ls {params}*_filtered.bed | grep '{wildcards.rbp}_' | while read line; "
        "do "
            "cat ${{line}} >> {output}; "
        "done"

rule sort_rbp: # sort RBP data
    input:
        rules.merge_cell_rbp.output
    output:
        temp(os.path.join(config["temp"],"{rbp}_sort.bed"))
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule rbp_formatted_dir: # create directory to store final RBP data & log files
    output:
        rbp = config["data"]["rbp_formatted"],
        log = config["logs"]
    shell:
        "mkdir -p {output.rbp} && mkdir -p {output.log}"

rule bedtools_merge_rbp: # merge concatenated & sorted RBP data
    input:
        rules.sort_rbp.output
    output:
        temp(os.path.join(config["temp"],"{rbp}_bedtools_merge.bed"))
    params:
        extra="-s -c 4,5,6,7,8 -o distinct,min,distinct,min,min"
    log:
        os.path.join(config["logs"],"bedtools_merge_rbp","{rbp}.log")
    wrapper:
        "v1.3.2/bio/bedtools/merge"

rule rbp_final_sort: # remove redundancy in bed file with gene names
    input:
        rules.bedtools_merge_rbp.output,
    output:
        os.path.join(config["data"]["rbp_formatted"],"{rbp}_uniq_regions.bed")
    shell:
        "cut -f1-6 {input} | sort -k1,1 -k2,3n -u > {output}"