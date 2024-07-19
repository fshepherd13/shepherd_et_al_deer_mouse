rule adapter_quality_trim:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1 = "../results/trimmed_reads/{id}_{rep}{pool}_R1_trimmed.fastq",
        r2 = "../results/trimmed_reads/{id}_{rep}{pool}_R2_trimmed.fastq"
    params:
        adapters = config["adapters"]
    log:
        "logs/trimming/{id}_{rep}{pool}.log"
    conda:
        "../envs/bbduk.yaml"
    shell:
        """
        bbduk.sh -Xmx1g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} minlen=50 qtrim=rl trimq=30 ktrim=r k=25 mink=11 ref="{params.adapters}" hdist=1 tpe tbo &> {log}
        """