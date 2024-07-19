rule star:
    input:
        fq1 = rules.trim.output.r1_paired,
        fq2 = rules.trim.output.r2_paired
    params:
        index = config["reference_genome"],
        outdir = '../results/star/{sample}'
    output:
        un_fq1 = '../results/star/{sample}/{sample}_Unmapped.out.mate1',
        un_fq2 = '../results/star/{sample}/{sample}_Unmapped.out.mate2',
    log:
        "logs/star/{sample}.log"
    conda:
        "../envs/star.yml"
    threads: 16
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.outdir}/{wildcards.sample}_ &> {log}
            
        rm -rf {params.outdir}/*_STARtmp
        """