rule merge_snoglobe_htrri:
    message: "Merge snoGloBe prediction and HTRRI interactions."
    input:
        os.path.join(config["data"]["snoglobe"],"pred_{sno}.98_3.gene.tsv")
    output:
        temp(os.path.join(config["filtered_data"]["sno_formatted"],"{sno}_snoglobe_htrri.bed"))
    params:
        config["data"]["HTRRI"],
        config["data"]["snoDB"],
        "{sno}"
    script:
        "../scripts/merge_snoglobe_htrri.py"

rule sort_sno:
    message: "Sort sno files for bedtools merge."
    input:
        rules.merge_snoglobe_htrri.output
    output:
        temp(os.path.join(config["filtered_data"]["sno_formatted"],"{sno}_sorted.bed"))
    shell:
        "awk \'{{print $1\"\t\"$2\"\t\"$3\"\t\"$12\"\t\"$5\"\t\"$6}}\' {input} | sort -k1,1 -k2,2n > {output}"

rule bedtools_merge_sno:
    message: "bedtools merge overlapping interactions."
    input:
        rules.sort_sno.output
    output:
        os.path.join(config["filtered_data"]["sno_formatted"],"{sno}.bed")
    params:
        extra="-s -c 4,5,6 -o distinct,mean,distinct"
    log:
        os.path.join(config["logs"],"bedtools_merge_sno","{sno}.log")
    wrapper:
        "v1.3.2/bio/bedtools/merge"