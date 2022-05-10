__author__ = "Kristina Song"
__email__ = "kristina.song@usherbrooke.ca"

import os
import glob
import pandas as pd

configfile:"config.json" # store paths and user-specific variables


#### global wildcard ####

ENCODE, = glob_wildcards(os.path.join(config["data"]["rbp_ENCODE"],"{encode}","replicate1.bed.gz"))



#### load snoRNA and RBP list ####

sno_list = pd.read_table(config["data"]["snoRNA_list"]).id.values.tolist()
rbp_list = pd.read_table(config["data"]["rbp_list"]).name1.values.tolist()



#### load rules ####

include: "rules/format_rbp.smk"
include: "rules/format_snoglobe.smk"
include: "rules/htrri.smk"
include: "rules/snodb.smk"


rule all:
    input:
        expand(os.path.join(config["data"]["rbp_formatted"],"{rbp}_uniq_regions.bed"),rbp=rbp_list),
        #os.path.join(config["outpath"],"snoglobe_targets.tsv"),
        expand(os.path.join(config["data"]["snoglobe_formatted"],"{sno}_uniq_regions.bed"), sno=sno_list),
        os.path.join(config["outpath"],"filtered_HTRRI.tsv"),
        os.path.join(config["outpath"],"snoDB_rbp_as_host_gene.tsv")
        #os.path.join(config["outpath"],"interaction_count.tsv")
"""
rule merge_interaction_count_files:
    # Merge all interaction counts into one file 
    input:
    output:
        os.path.join(config["outpath"],"interaction_count.tsv")
    shell:
        "touch {output}"

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
"""