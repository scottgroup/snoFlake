__author__ = "Kristina Song"
__email__ = "kristina.song@usherbrooke.ca"

import os
import glob

configfile:"config.json" # store paths and variables

ENCODE, = glob_wildcards(os.path.join(config["data"]["rbp_ENCODE"],"{encode}","replicate1.bed.gz"))

rule all:
    input:
        expand(os.path.join(config["data"]["rbp_formatted"],"{rbp}_uniq_regions.bed"),rbp=config["rbp"]),
        os.path.join(config["outpath"],"sno_embedded_in_rbp_host_gene.tsv"),
        os.path.join(config["outpath"],"snoglobe_targets.tsv"),
        expand(os.path.join(config["data"]["snoglobe_formatted"],"{sno}_unique_regions.bed"), sno=config["snoRNA"])

rule filter_rbp: # filter raw ENCODE RBP data by p-value -log10(p-val) >= 5 and by fold-enrichment log2(fold-enrichment) >= 3
    input:
        os.path.join(config["data"]["rbp_ENCODE"],"{encode}","replicate{rep}.bed.gz")
    output:
        temp(os.path.join(config["data"]["temp"],"{encode}_replicate{rep}_filtered.bed"))
    params:
        temp_dir = config["data"]["temp"],
        p_thres = config["p_threshold"] # specify p-value threshold in config.json
    shell:
        "mkdir -p {params.temp_dir} && "
        "gunzip -c {input}| sort -k8 -n | awk -F'\t' '$8>={params.p_thres} && $7>=3' > {output}"

rule merge_rep_rbp: # merge RBP replicate data
    input:
        rep1 = os.path.join(config["data"]["temp"],"{encode}_replicate1_filtered.bed"),
        rep2 = os.path.join(config["data"]["temp"],"{encode}_replicate2_filtered.bed")
    output:
        temp(os.path.join(config["data"]["temp"],"{encode}_rep_merge.bed"))
    shell:
        "cat {input.rep1} >> {output} && cat {input.rep2} >> {output}"

rule create_empty_merge_file: # empty temp merge files to automatically delete them after run
    output:
        temp(os.path.join(config["data"]["temp"],"{rbp}_all_merge_empty.bed"))
    shell:
        "touch {output}"

rule merge_cell_rbp_temp: # merge by cell line
    input:
        empty = rules.create_empty_merge_file.output,
        dir = config["data"]["rbp_ENCODE"],
        encode = expand(rules.merge_rep_rbp.output, encode=ENCODE)
    output:
        temp(os.path.join(config["data"]["temp"],"{rbp}_all_merge_temp.bed"))
    params:
        grep_rbp = "{rbp}_",
        temp_dir = config["data"]["temp"]
    shell:
        "ls {input.dir} | grep '{params.grep_rbp}' | sed 's/{params.grep_rbp}//g' | while read line;"
        "do "
            "cat {params.temp_dir}{wildcards.rbp}_${{line}}_rep_merge.bed >> {output};"
        "done"

rule merge_cell_rbp: # move merged RBP data from temp to final file
    input:
        rules.merge_cell_rbp_temp.output
    output:
       temp(os.path.join(config["data"]["temp"],"{rbp}_all_merge.bed"))
    shell:
        "mv {input} {output}"

rule sort_rbp: # sort RBP data
    input:
        rules.merge_cell_rbp.output
    output:
        temp(os.path.join(config["data"]["temp"],"{rbp}_sort.bed"))
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule rbp_formatted_dir: # create directory to store final RBP data & logs
    output:
        rbp = config["data"]["rbp_formatted"],
        log = config["data"]["logs"]
    shell:
        "mkdir -p {output.rbp} && mkdir -p {output.log}"

rule bedtools_merge_rbp: # merge concatenated & sorted RBP data
    input:
        rules.sort_rbp.output
    output:
        temp(os.path.join(config["data"]["temp"],"{rbp}_bedtools_merge.bed"))
    params:
        extra="-s -c 4,5,6,7,8 -o distinct,min,distinct,min,min"
    log:
        os.path.join(config["data"]["logs"],"{rbp}_merge.log")
    wrapper:
        "v1.3.1/bio/bedtools/merge"

rule rbp_final_sort: # remove redundancy in bed file with gene names
    input:
        rules.bedtools_merge_rbp.output
    output:
        os.path.join(config["data"]["rbp_formatted"],"{rbp}_uniq_regions.bed")
    shell:
        "cut -f1-6 {input} | sort -k1,1 -k2,3n -u > {output}"

rule sno_embedded_in_host_gene: # find snoRNAs that have protein coding host genes and keep those that belong to our RBP group
    input:
        snorna = config["data"]["snoRNA_list"],
        rbp = config["data"]["rbp_list"],
        snodb = config["data"]["snoDB_host_gene"]
    output:
        os.path.join(config["outpath"],"sno_embedded_in_rbp_host_gene.tsv")
    shell:
        "echo \"snoRNA\tRBP\tinteraction\"> {output} && "
        "cat {input.snorna} | while read line; do if grep -q $line {input.snodb}; then "
        "result=$(grep $line {input.snodb}); rbp=$(awk '{{print $2}}' <<< $result); "
        "if grep -q $rbp {input.rbp}; then echo $line\"\t\"$rbp\"\t\"\"embedded_sno\" >> {output}; fi; fi; done"

rule snoglobe_targets: # get all unique targets for each snoRNA from snoGloBe predictions
    input:
        config["data"]["snoRNA_list"]
    output:
        os.path.join(config["outpath"],"snoglobe_targets.tsv")
    params:
        config["data"]["snoglobe"]
    shell:
        "echo \"snoRNA\ttarget\tinteraction\"> {output} && "
        "cat {input} | while read line; do cat {params}pred_$line.98_3.gene.tsv | cut -f 4,12 | sort -u | awk '{{print $1\"\t\"$2\"\t\"\"snoglobe\"}}' >> {output}; done "

rule snoglobe_uniq: # format snoglobe prediction for snoRNA-RBP p-val calculation
    input:
        os.path.join(config["data"]["snoglobe"],"pred_{sno}.98_3.gene.tsv")
    output:
        os.path.join(config["data"]["snoglobe_formatted"],"{sno}_unique_regions.bed")
    params:
        config["data"]["snoglobe_formatted"]
    shell:
        "mkdir -p {params} && cut -f1-6 {input} | sort -k1,1 -k2,3n -u > {output}"
