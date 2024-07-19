rule concatenate:
    input:
        r1 = expand("../results/star/{sample}/{sample}_Unmapped.out.mate1", sample=SAMPLES),
        r2 = expand("../results/star/{sample}/{sample}_Unmapped.out.mate2", sample=SAMPLES)
    output:
        r1 = temp(expand("../results/star/{exp}_all_unmapped_1.fq", exp=EXP)),
        r2 = temp(expand("../results/star/{exp}_all_unmapped_2.fq", exp=EXP))
    shell:
        '''
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        '''
        
rule modify_for_trinity:
    input:
        r1 = expand("../results/star/{exp}_all_unmapped_1.fq", exp=EXP),
        r2 = expand("../results/star/{exp}_all_unmapped_2.fq", exp=EXP)
    output:
        r1 = temp(expand("../results/star/{exp}_all_unmapped_ed1.fq", exp=EXP)),
        r2 = temp(expand("../results/star/{exp}_all_unmapped_ed2.fq", exp=EXP))
    shell:
        '''
        #Append read direction information using awk
        #On every 4th line of the unmapped_1.fq file (NR%4==1), strip white space and add
        # "/1" text at the end. Print the lines to a new fq file.
        awk '{{if (NR%4==1) $1=$1 "/1"; print}}' {input.r1} > {output.r1}
        awk '{{if (NR%4==1) $1=$1 "/2"; print}}' {input.r2} > {output.r2}
        '''