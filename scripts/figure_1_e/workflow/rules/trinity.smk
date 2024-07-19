rule trinity:
    input:
        r1 = expand("../results/star/{exp}_all_unmapped_ed1.fq", exp=EXP),
        r2 = expand("../results/star/{exp}_all_unmapped_ed2.fq", exp=EXP)
    output:
        expand("../results/trinity/{exp}_unmapped_denovo_contigs.fasta", exp=EXP)
    params:
        scratch=expand("/scratch.global/sheph085/trinity_{exp}", exp=EXP)
    conda:
        "../envs/trinity.yml"
    log:
        "logs/trinity/Trinity.log"
    threads: 64
    shell:
        '''
        Trinity --seqType fq \
                --max_memory 480G \
                --left {input.r1} \
                --right {input.r2} \
                --SS_lib_type RF \
                --CPU {threads} \
                --min_contig_length 150 \
                --normalize_reads \
                --full_cleanup \
                --output {params.scratch} &> {log}
                
        cp {params.scratch}.Trinity.fasta {output}
        '''
