rule consensus_depth:
    input:
        "../results/consensus_calling/final/{pet_store}.merged.bam"
    output:
        "../results/depth_calc/consensus/{pet_store}.depth"
    params:
        index = config["ref_seq"]+".fa"
    conda:
        "../envs/ivar.yaml"
    shell:
        '''
        samtools depth {input} > {output}
        '''

rule depth_by_rep:
    input:
        "../results/variant_calling/merged_pools/{id}_{rep}.merged.bam"
    output:
        "../results/depth_calc/by_rep/{id}_{rep}.depth"
    params:
        index = expand("../results/consensus_calling/{pet_store}.fa", pet_store=PET_STORE)    
    conda:
        "../envs/ivar.yaml"
    shell:
        '''
        samtools depth {input} > {output}
        '''

rule depth_by_pool:
    input:
        "../results/variant_calling/aligned/{id}_{rep}{pool}.sorted.bam"
    output:
        "../results/depth_calc/by_pool/{id}_{rep}{pool}.depth"
    params:
        index = expand("../results/consensus_calling/{pet_store}.fa", pet_store=PET_STORE)
    conda:
        "../envs/ivar.yaml"
    shell:
        '''
        samtools depth {input} > {output}
        '''