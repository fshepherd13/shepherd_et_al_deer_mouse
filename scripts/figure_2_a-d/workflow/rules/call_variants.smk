#Align reads to the pet store consensus sequence
rule align_to_consensus:
    input:
        r1 = "../results/trimmed_reads/{id}_{rep}{pool}_R1_trimmed.fastq",
        r2 = "../results/trimmed_reads/{id}_{rep}{pool}_R2_trimmed.fastq"
    output:
        "../results/variant_calling_{pet_store}/aligned/{id}_{rep}{pool}.sorted.bam"
    params:
        index = "../results/consensus_calling/final/{pet_store}"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        mkdir -p $(dirname {output})/
        bwa mem {params.index} {input.r1} {input.r2} | samtools view -F 4 -Sb | samtools sort -T {wildcards.id}_{wildcards.rep}{wildcards.pool}_align -o {output}
        samtools index {output}
        """

rule merge_sample_pools:
    input:
        expand("../results/variant_calling_{{pet_store}}/aligned/{{id}}_{{rep}}{pool}.sorted.bam", pool=["1","2"])
    output:
        "../results/variant_calling_{pet_store}/merged_pools/{id}_{rep}.merged.bam"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        samtools merge {output} {input} 
        samtools index {output}
        """

rule realigned_trim_primer_quality:
    input:
        bam = "../results/variant_calling_{pet_store}/merged_pools/{id}_{rep}.merged.bam",
        bed = "../results/bed/{pet_store}.bed",
    output:
        "../results/variant_calling_{pet_store}/trimmed/{id}_{rep}.trimmed.sorted.bam"
    params:
        prefix = "../results/variant_calling_{pet_store}/trimmed/{id}_{rep}.trimmed"
    log:
        "logs/variant_calling_{pet_store}/{id}_{rep}.primer_trim.log"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        mkdir -p $(dirname {output})/
        ivar trim -b {input.bed} -p {params.prefix} -i {input.bam}
        samtools sort -T {wildcards.id}_{wildcards.rep}_trim -o {output} {params.prefix}.bam
        samtools index {output}
        """

rule remove_reads: #Get rid of reads associated with primers with iSNVs in binding site
    input:
        "../results/variant_calling_{pet_store}/trimmed/{id}_{rep}.trimmed.sorted.bam",
        "../results/{pet_store}_masked_primer_names.txt",
        "../results/bed/{pet_store}.bed"
    output:
        "../results/variant_calling_{pet_store}/masked/{id}_{rep}.masked.sorted.bam"
    params:
        prefix = "../results/variant_calling_{pet_store}/masked/{id}_{rep}.masked"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        mkdir -p $(dirname {output})/
        ivar removereads -i {input[0]} -p {params.prefix} -t {input[1]} -b {input[2]}
        samtools sort -T ../results/variant_calling/{wildcards.id}_{wildcards.rep}_remove_reads -o {output} {params.prefix}.bam
        samtools index {output}
        """

rule call_variants_post_removal:
    input:
        bam = "../results/variant_calling_{pet_store}/masked/{id}_{rep}.masked.sorted.bam",
        ref = "../results/consensus_calling/final/{pet_store}.fa",
        gff = "../results/consensus_calling/final/{pet_store}.gff"
    output:
        "../results/variant_calling_{pet_store}/unfiltered/{id}_{rep}.tsv"
    conda:
        "../envs/ivar.yaml"
    log:
        "logs/variant_calling_{pet_store}/unfiltered.{id}_{rep}.log"
    shell:
        """
        samtools mpileup -A -d 0 --reference {input.ref} -Q 0 {input.bam} | ivar variants -p {output} -m 100 -q 30 -t 0.03 -g {input.gff} -r {input.ref} &> {log}
        """

#Compare the variants that are called between each of the replicates and filter the ones that only occur in one replicate to avoid false positive variant calling.
rule filter_variants:
    input:
        expand("../results/variant_calling_{{pet_store}}/unfiltered/{{id}}_{rep}.tsv", rep=["A","B"])
    output:
        "../results/variant_calling_{pet_store}/final/{id}.filtered.tsv"
    conda:
        "../envs/ivar.yaml"
    log:
        "logs/variant_calling_{pet_store}/final/filtered.{id}.log",
    shell:
        """
        mkdir -p $(dirname {output})/
        ivar filtervariants -p {output} {input} &> {log}
        """