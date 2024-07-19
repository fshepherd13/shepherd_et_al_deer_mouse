rule align_petstore_reads:
    input:
        r1 = "../results/trimmed_reads/{pet_store}_{rep}{pool}_R1_trimmed.fastq",
        r2 = "../results/trimmed_reads/{pet_store}_{rep}{pool}_R2_trimmed.fastq",
        ref = expand("../ref_files/"+config["ref_seq"]+".{ext}", ext=["amb","ann","bwt","pac","sa"])
    params:
        index=config["ref_seq"]
    log:
        "logs/generate_consensus/{pet_store}_{rep}{pool}.align.log"
    output:
        "../results/consensus_calling/aligned/{pet_store}_{rep}{pool}.sorted.bam"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        mkdir -p $(dirname {output})/
        bwa mem {params.index} {input.r1} {input.r2} | samtools view -F 4 -b | samtools sort -T {wildcards.pet_store}_align -o {output} &> {log}
        samtools index {output}
        """

rule merge_pools:
    input:
        expand("../results/consensus_calling/aligned/{{pet_store}}_{{rep}}{pool}.sorted.bam", pool = ["1", "2"])
    log:
        "logs/merge_pools/{pet_store}_{rep}.align.log"
    output:
        "../results/consensus_calling/merged_pools/{pet_store}_{rep}.sorted.bam"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        samtools merge {output} {input}
        samtools index {output}
        """

rule trim_primer_quality:
    input:
        bam = "../results/consensus_calling/merged_pools/{pet_store}_{rep}.sorted.bam",
        bed = config["primers"]
    params:
        prefix = "../results/consensus_calling/trimmed/{pet_store}_{rep}.trimmed"
    log:
        "logs/consensus_calling/{pet_store}_{rep}.primer_trim.log"
    conda:
        "../envs/ivar.yaml"
    output:
        "../results/consensus_calling/trimmed/{pet_store}_{rep}.trimmed.sorted.bam"
    shell:
        """
        mkdir -p $(dirname {output})/
        ivar trim -b {input.bed} -p {params.prefix} -i {input.bam}
        samtools sort -T {wildcards.pet_store}_{wildcards.rep}_trim -o {output} {params.prefix}.bam
        samtools index {output}
        """

rule merge_petstore_replicates:
    input:
        expand("../results/consensus_calling/trimmed/{{pet_store}}_{rep}.trimmed.sorted.bam", rep = ["A", "B"])
    output:
        "../results/consensus_calling/final/{pet_store}.merged.bam"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        samtools merge {output} {input}
        samtools index {output}
        """

rule call_petstore_consensus:
    input:
        "../results/consensus_calling/final/{pet_store}.merged.bam"
    output:
        fa = "../results/consensus_calling/final/{pet_store}.fa",
        index = expand("../results/consensus_calling/final/{{pet_store}}"+".{ext}", ext=["amb","ann","bwt","pac","sa"])
    log:
        pileup = "logs/generate_consensus/{pet_store}.call_consensus.log",
        index = "logs/generate_consensus/{pet_store}.index_consensus.log"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        samtools mpileup -aa -A -d 0 -Q 0 {input} | ivar consensus -q 30 -p {output.fa} -m 10 &> {log.pileup}
        bwa index -p ../results/consensus_calling/final/{wildcards.pet_store} {output.fa} &> {log.index}
        """

