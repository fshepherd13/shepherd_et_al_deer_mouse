rule create_primer_fasta:
    input:
        "{csv}".format(csv = config["primer_csv"])
    output:
        "../results/primers/primers.fa"
    shell:
        "cut -f 1,2 -d ',' {input} | sed 's/^/>/g' | tr ',' '\n' > {output}"
        
rule create_primer_bam: #Aligns primers to the pet store consensus sequence, outputs bam file. bwa mem allows for some mismatches between primer and reference seq but not a lot, so not all primers will align. 
    input:
        index = expand("../results/consensus_calling/final/{{pet_store}}"+".{ext}", ext=["amb","ann","bwt","pac","sa"]),
        fa = "../results/primers/primers.fa"
    output:
        "../results/bed/{pet_store}.bam"
    params:
        index= "../results/consensus_calling/final/{pet_store}"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        bwa mem -k 5 -T 16 -B 0.5 {params.index} {input.fa} | samtools view -b -F 4 | samtools sort -T ../results/bed/{wildcards.pet_store}_align -o {output}
        samtools index {output}
        """
        
rule create_bed: #Creates bed file of positions where primers align to the pet store consensus sequence
    input:
        "../results/bed/{pet_store}.bam"
    output:
        "../results/bed/{pet_store}.bed"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        bedtools bamtobed -i {input} > {output}
        """
        
rule call_variants_in_primer: #Find places where there are iSNVs called in primer binding sites against the generated pet store consensus sequence. These primers will be masked and associated reads removed from variant calling
    input:
        bam = "../results/bed/{pet_store}.bam",
        fa = "../results/consensus_calling/final/{pet_store}.fa"
    output:
        "../results/primers/{pet_store}_primer_mismatches.tsv"
    conda:
        "../envs/ivar.yaml"
    shell:
        """
        samtools mpileup -A -B -d 0 --reference {input.fa} -Q 0 {input.bam} | ivar variants -p {output} -t 0.03
        """

rule get_masked: #
    input:
        "../results/primers/{pet_store}_primer_mismatches.tsv",
        "../results/bed/{pet_store}.bed",
        "{pair_information}".format(pair_information = config["pair_information"])
    output:
        "../results/{pet_store}_masked_primer_names.txt"
    conda:
        "../envs/ivar.yaml"
    shell:
        "ivar getmasked -i {input[0]} -b {input[1]}  -f {input[2]} -p {output}"