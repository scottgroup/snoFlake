__author__ = "Kristina Song"
__email__ = "kristina.song@usherbrooke.ca"

import os
import glob
import pandas as pd

configfile:"config.json" # store paths and user-specific variables


#### global wildcard ####

#ENCODE, = glob_wildcards(os.path.join(config["data"]["rbp_ENCODE"],"{encode}","replicate1.bed.gz"))



#### load snoRNA and RBP list ####
# Need to get "nodes" output before importing these files
#sno_list = pd.read_table(config["data"]["snoRNA_list"]).id.values.tolist()
#rbp_list = pd.read_table(config["data"]["rbp_list"]).RBP1.values.tolist()



#### load rules ####

include: "rules/rank_sno_rbp.smk"
#include: "rules/format_rbp.smk"
#include: "rules/format_snoglobe.smk"
#include: "rules/htrri.smk"
#include: "rules/snodb.smk"
#include: "rules/binding.smk"
#include: "rules/string.smk"
#include: "rules/overlaps.smk"


rule all:
    input:
        # Filtering and sorting snoRNAs and RBPs for node selection
        # Run FIRST!!!
        nodes = {
            config["data"]["snoRNA_list"],
            config["data"]["rbp_list"],
            os.path.join(config["outpath"],"snoRNA_ranking.tsv"),
            os.path.join(config["outpath"],"rbp_ranking.tsv")
        }
        #interactions = { # interactions to be included in network
        #    os.path.join(config["outpath"],"rbp_bind_to_sno_transcript.tsv")
        #}
        #expand(os.path.join(config["data"]["rbp_formatted"],"{rbp}_uniq_regions.bed"),rbp=rbp_list),
        #expand(os.path.join(config["data"]["snoglobe_formatted"],"{sno}_uniq_regions.bed"), sno=sno_list),
        #os.path.join(config["outpath"],"filtered_HTRRI.tsv"),
        #os.path.join(config["outpath"],"snoDB_rbp_as_host_gene.tsv"),
        #os.path.join(config["outpath"],"snoglobe_sno_bind_to_rbp_transcript.tsv"),
        #os.path.join(config["outpath"],"STRING_physical_binding.tsv"),
        #os.path.join(config["outpath"],"significant_snoglobe_sno_rbp_target_overlaps.tsv")
        #os.path.join(config["outpath"],"interaction_counts.tsv"),
        #expand(os.path.join(config["outpath"],"sno_rbp_overlaps_p_vals","{sno}_rbp_overlaps.tsv"),sno=sno_list),
        #expand(os.path.join(config["outpath"],"sno_sno_overlaps_p_vals","{sno}_sno_overlaps.tsv"),sno=sno_list),
        #expand(os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals","{rbp}_rbp_overlaps.tsv"),rbp=rbp_list)
        #os.path.join(config["outpath"],"rbp_rbp_overlaps_p_vals","NSUN2_rbp_overlaps.tsv")
        #os.path.join(config["outpath"],"significant_snoglobe_sno_rbp_target_overlaps.tsv"),
        #os.path.join(config["outpath"],"significant_rbp_rbp_target_overlaps.tsv"),
        #os.path.join(config["outpath"],"significant_sno_sno_target_overlaps.tsv")
"""
rule merge_interaction_count_files:
    message: "Merge all interaction counts into one file."
    input:
        string = rules.merge_all_STRING_counts.output,
        host_gene = rules.filter_snodb_host_gene.output.counts,
        htrri = rules.get_htrri.output.counts,
        sno_rbp_transcript = rules.sno_rbp_transcript_out_file.output.counts,
        rbp_sno_transcript = rules.rbp_sno_transcript_out_file.output.counts
    output:
        os.path.join(config["outpath"],"interaction_counts.tsv")
    shell:
        "echo -e \"INTERACTION_TYPE\t$(basename {input.string})\" >> {output} && cat {input.string} >> {output}; "
        "echo -e \"INTERACTION_TYPE\t$(basename {input.host_gene})\" >> {output} && cat {input.host_gene} >> {output}; "
        "echo -e \"INTERACTION_TYPE\t$(basename {input.htrri})\" >> {output} && cat {input.htrri} >> {output}; "
        "echo -e \"INTERACTION_TYPE\t$(basename {input.sno_rbp_transcript})\" >> {output} && cat {input.sno_rbp_transcript} >> {output}; "
        "echo -e \"INTERACTION_TYPE\t$(basename {input.rbp_sno_transcript})\" >> {output} && cat {input.rbp_sno_transcript} >> {output}; "
"""