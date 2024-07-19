rule fastqc_raw:
    input:
        r1=get_r1,
        r2=get_r2
    output:
        r1="qc/untrimmed/{sample}_R1_fastqc.html",
        r2="qc/untrimmed/{sample}_R2_fastqc.html"
    shell:
        """
        #!/bin/bash
        fastqc {input.r1} {input.r2} -q -o qc/
        """

rule fastqc_trimmed:
    input:
        r1="../results/trimmed/{sample}_R1_trimmed_paired.fq.gz",
        r2="../results/trimmed/{sample}_R2_trimmed_paired.fq.gz"
    output:
        r1="qc/untrimmed/{sample}_R1_fastqc.html",
        r2="qc/untrimmed/{sample}_R2_fastqc.html"
    shell:
        """
        #!/bin/bash
        fastqc {input.r1} {input.r2} -q -o qc/
        """