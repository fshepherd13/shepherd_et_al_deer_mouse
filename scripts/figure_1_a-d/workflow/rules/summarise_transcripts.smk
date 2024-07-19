rule summarize_transcripts:
    input:
        lineages = rules.assign_lineages.output,
        salmon_files = expand("../results/salmon_quant/{sample}/quant.sf", sample=SAMPLES)
    output:
        species_counts = expand("../results/final/{exp}_sp_counts.csv", exp=EXP),
        family_counts = expand("../results/final/{exp}_family_counts.csv", exp=EXP)
    log:
        "logs/summarise_transcripts.log"
    conda:
        "../envs/deseq2.yml"
    script:
        "../scripts/deseq2.R"