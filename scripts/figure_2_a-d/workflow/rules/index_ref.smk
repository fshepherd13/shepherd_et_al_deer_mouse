rule index_ref:
    input:
        config["ref_seq"]+".fa"
    output:
        expand("../ref_files/"+config["ref_seq"]+".{ext}", ext=["amb","ann","bwt","pac","sa"])
    params:
        prefix=config["ref_seq"]
    conda:
        "../envs/ivar.yaml"
    shell:
        '''
        bwa index {input} -p {params.prefix}
        '''
