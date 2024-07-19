rule call_snvs:
    input:
        bam = config["in_dir"]+"{id}_{rep}.merged.bam"
    output:
        "../results/snvs/{id}_{rep}.vcf.gz"
    params:
        ref = config["ref_seq"]+".fa"
    log:
        "logs/call_snvs/{id}_{rep}.log"
    shell:
        '''
        module load bcftools
        bcftools mpileup -A -f {params.ref} -d 100000000 --annotate FORMAT/AD,FORMAT/DP {input.bam} | bcftools call -mv -Oz -o {output} &> {log}
        '''