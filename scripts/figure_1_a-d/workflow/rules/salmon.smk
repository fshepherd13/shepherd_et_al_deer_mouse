rule salmon_index:
    input:
        rules.trinity.output
    output:
        directory("../results/salmon_indices/salmon_unmapped_index")
    log:
        "logs/salmon/index.log"
    threads: 16
    conda:
        "../envs/salmon.yml"
    shell:
        """
        salmon index -t {input} -i {output} -p 16 &> {log}
        """
        
rule salmon_quant:
    input:
        r1 = rules.trim.output.r1_paired,
        r2 = rules.trim.output.r2_paired,
        index = "../results/salmon_indices/salmon_unmapped_index"
    output:
        "../results/salmon_quant/{sample}/quant.sf"
    params:
        dir = directory("../results/salmon_quant/{sample}/")
    threads: 16
    conda:
        "../envs/salmon.yml"
    log:
        "logs/salmon/{sample}.quant.log"
    shell:
        """
        salmon quant -i {input.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            --seqBias \
            --gcBias \
            -o {params.dir} &> {log}
        """