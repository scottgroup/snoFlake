rule filter_merge_ENCODE:
    output:
        "results/interactions/ENCODE/{rbp}_filtered_merged.bed"
    params:
        in_dir = config["path"]["ENCODE"],
        tmp_file = "results/interactions/ENCODE/{rbp}_filtered_merged_tmp.bed",
        thres = config["thresholds"]["ENCODE"]
    message:
        "Filter raw ENCODE {wildcards.rbp} data by p-value and signal value thresholds of {params.thres} and merge replicates and cell lines."
    shell:
        "ls {params.in_dir} | grep '{wildcards.rbp}_' | while read line; "
        "do "
            "gunzip -c {params.in_dir}/${{line}}/replicate1.bed.gz | sort -k8 -n | awk -F'\t' '$8>={params.thres} && $7>={params.thres}' >> {params.tmp_file}; "
            "gunzip -c {params.in_dir}/${{line}}/replicate2.bed.gz | sort -k8 -n | awk -F'\t' '$8>={params.thres} && $7>={params.thres}' >> {params.tmp_file}; "
        "done && "
        "awk \'{{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$8\"\t\"$6}}\' {params.tmp_file} | sort -k1,1 -k2,2n > {output} && "
        "rm -f {params.tmp_file} && "
        "sed -i \'/random/d\' {output} && "
        "sed -i \'/Un_/d\' {output}"
        

rule bedtools_merge_ENCODE:
    input:
        rules.filter_merge_ENCODE.output
    output:
        "results/interactions/ENCODE/{rbp}.bed"
    conda:
        "../envs/bedtools.yaml"
    message:
        "Bedtools merge ENCODE {wildcards.rbp} interactions."
    shell:
        "bedtools merge -s -c 4,5,6 -o distinct,min,distinct -i {input} | sort -k1,1 -k2,2n -u | awk -v var=\"{wildcards.rbp}\" \'{{print $1\"\t\"$2\"\t\"$3\"\t\"var\"\t\"$5\"\t\"$6}}\' | sed 's/chrM/chrMT/g' | sed 's/^chr\|%$//g' > {output} && "
        "rm -f {input}"


rule bedtools_merge_snoGloBe_HTRRI:
    input:
        snoglobe = os.path.join(config["path"]["snoglobe"],"pred_{sno}.95_3.gene.tsv"),
	    htrri = config["path"]["HTRRI"]
    output:
        "results/interactions/snoGloBe_HTRRI/{sno}.bed"
    params:
        snoglobe_thres = config["thresholds"]["snoGloBe"],
        support_thres = config["thresholds"]["HTRRI"]
    conda:
        "../envs/bedtools.yaml"
    message:
        "Format snoRNA interactions obtained from snoGloBe predictions (-t {params.snoglobe_thres}) and HTRRI (-t {params.support_thres}) for {wildcards.sno} to run bedtools merge."
    shell:
        "python3 workflow/scripts/merge_snoglobe_htrri.py {input.snoglobe} {input.htrri} {wildcards.sno} {output} {params.snoglobe_thres} {params.support_thres}"
