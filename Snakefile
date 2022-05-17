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
include: "rules/binding.smk"
include: "rules/string.smk"
#include: "rules/overlaps.smk"


rule all:
    input:
        expand(os.path.join(config["data"]["rbp_formatted"],"{rbp}_uniq_regions.bed"),rbp=rbp_list),
        expand(os.path.join(config["data"]["snoglobe_formatted"],"{sno}_uniq_regions.bed"), sno=sno_list),
        os.path.join(config["outpath"],"filtered_HTRRI.tsv"),
        os.path.join(config["outpath"],"snoDB_rbp_as_host_gene.tsv"),
        os.path.join(config["outpath"],"sno_bind_to_rbp_transcript.tsv"),
        os.path.join(config["outpath"],"rbp_bind_to_sno_transcript.tsv"),
        os.path.join(config["outpath"],"STRING_physical_binding.tsv"),
        os.path.join(config["outpath"],"interaction_counts.tsv")
        #expand(os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals","{sno}.tsv"),sno=sno_list),
        #expand(os.path.join(config["outpath"],"sno_sno_overlaps_p_vals","{sno}.tsv"),sno=sno_list),
        #expand(os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals","{rbp}.tsv"),rbp=rbp_list)

rule merge_interaction_count_files:
    message: "Merge all interaction counts into one file."
    input:
        string = rules.merge_all_STRING_counts.output,
        host_gene = rules.filter_snodb_host_gene.output.counts,
        htrri = rules.get_htrri.output.counts
    output:
        os.path.join(config["outpath"],"interaction_counts.tsv")
    shell:
        "echo -e \"INTERACTION_TYPE\t$(basename {input.string})\" >> {output} && cat {input.string} >> {output}; "
        "echo -e \"INTERACTION_TYPE\t$(basename {input.host_gene})\" >> {output} && cat {input.host_gene} >> {output}; "
        "echo -e \"INTERACTION_TYPE\t$(basename {input.htrri})\" >> {output} && cat {input.htrri} >> {output}; "