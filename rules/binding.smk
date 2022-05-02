rule sno_rbp_transcript: # snoRNA binds to RBP transcript
    input:
    output:
    parmas:
        
    log:
        os.path.join(config["logs"],"sno_rbp_transcript","{sno}_intersect.log")
    wrapper:
        "v1.3.2/bio/bedtools/intersect"





    input:
        rules.sort_rbp.output
    output:
        temp(os.path.join(config["temp"],"{rbp}_bedtools_merge.bed"))
    params:
        extra="-s -c 4,5,6,7,8 -o distinct,min,distinct,min,min"
    log:
        os.path.join(config["logs"],"{rbp}_merge.log")
    wrapper:
        "v1.3.1/bio/bedtools/merge"