rule fastqc_trimmed:
    input:
        rules.adapter_quality_trim.output.r1,
        rules.adapter_quality_trim.output.r2
    output:
        "../results/quality/{sample}_R1_trimmed_fastqc.html",
        "../results/quality/{sample}_R2_trimmed_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    shell:
        '''
        #!/bin/bash
        fastqc {input} -q -o ..results/quality/
        '''