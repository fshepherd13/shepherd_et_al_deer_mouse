rule fa_split:
    input:
        rules.trinity.output
    output:
        temp(expand("../results/trinity/{exp}_unmapped_denovo_contigs.{num}.fasta", exp=EXP, num=NUM)),
        temp(expand("../results/trinity/{exp}_unmapped_denovo_contigs.fasta.flat", exp=EXP)),
        temp(expand("../results/trinity/{exp}_unmapped_denovo_contigs.fasta.gdx", exp=EXP))
    conda:
        "../envs/blast.yml"
    shell:
        "pyfasta split {input} -n 100"

rule blast:
    input:
        expand("../results/trinity/{exp}_unmapped_denovo_contigs.{{num}}.fasta", exp=EXP)
    output:
        temp(expand("../results/blast_results/{exp}.{{num}}.out", exp=EXP))
    params:
        db = config["blast"]["db"]
    conda:
        "../envs/blast.yml"
    threads: 8
    log:
        "logs/blast/{{num}}_blast.log"
    shell:
        '''
        blastn -query {input} -db {params.db} -out {output} -max_target_seqs 10 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -num_threads 8 &> {log}
        '''
