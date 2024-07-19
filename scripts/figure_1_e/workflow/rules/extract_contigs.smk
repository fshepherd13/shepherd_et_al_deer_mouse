rule extract_contigs:
    input:
        lineages = "../results/lineages/{sample}_transcript_lineages.csv",
        trinity = "../results/trinity/{sample}.Trinity.fasta"
    output:
        contig_ids = expand("../results/viral_contigs/{{sample}}/{{sample}}_{virus_family}_contigs.txt", virus_family=VIRUS_FAMILY),
        contig_fasta = expand("../results/viral_contigs/{{sample}}/{{sample}}_{virus_family}_contigs.fasta", virus_family=VIRUS_FAMILY)
    params:
        virus_family=VIRUS_FAMILY,
        out_dir = "../results/viral_contigs/{sample}"
    shell:
        '''
        mkdir -p {params.out_dir}
        grep -i {params.virus_family} {input.lineages} > {output.contig_ids} || true
        cat {output.contig_ids} | cut -d ',' -f1 | grep -F -A 1 -f - {input.trinity} > {output.contig_fasta} || true
        '''
        
rule aggregate_contigs:
    input:
        expand("../results/viral_contigs/{sample}/{sample}_{virus_family}_contigs.fasta", sample = SAMPLES, virus_family = VIRUS_FAMILY)
    output:
        expand("../results/viral_contigs/{virus_family}_contigs.fasta", virus_family = VIRUS_FAMILY)
    shell:
        '''
        tail -n +1 {input} > {output}
        '''
    