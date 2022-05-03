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
rbp_list = pd.read_table(config["data"]["rbp_list"]).name.values.tolist()



#### load rules ####

include: "rules/format_rbp.smk"
include: "rules/format_snoglobe.smk"
include: "rules/htrri.smk"


rule all:
    input:
        expand(os.path.join(config["data"]["rbp_formatted"],"{rbp}_uniq_regions.bed"),rbp=rbp_list),
        #os.path.join(config["outpath"],"sno_embedded_in_rbp_host_gene.tsv"),
        #os.path.join(config["outpath"],"snoglobe_targets.tsv"),
        expand(os.path.join(config["data"]["snoglobe_formatted"],"{sno}_uniq_regions.bed"), sno=sno_list),
        os.path.join(config["outpath"],"filtered_HTRRI.tsv")

"""
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
"""