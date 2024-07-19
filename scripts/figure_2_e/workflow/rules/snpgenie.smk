rule reheader:
    input:
        vcf = "../results/filtered_snvs/{id}/{id}_filtered_merged.vcf"
    output:
        vcf_new = "../results/snpgenie/{id}_filtered_merged_ed.vcf"
    params:
        header_file = "../results/filtered_snvs/{id}/headers.txt",
        tempfile = "../results/filtered_snvs/{id}/{id}_filtered_merged_ed.vcf",
        dir = "../results/snpgenie/"
    conda:
        "../envs/env.yaml"
    message: "Editing header of vcf file so that it works with SNPGenie (removes the .merged.bam suffix)"
    shell:
        '''
        module load bcftools
        echo {wildcards.id}_A > {params.header_file}
        echo {wildcards.id}_B >> {params.header_file}
        bcftools reheader -s {params.header_file} {input.vcf} -o {params.tempfile}
        mv {params.tempfile} {params.dir}
        touch {output.vcf_new}
        '''


rule snpgenie:
    input:
        ref = config["ref_seq"]+".fa",
        gtf = config["ref_seq"]+".gtf",
        vcf = rules.reheader.output.vcf_new
    params:
        dir = "/home/langlois/sheph085/dirty_mouse_virome/fsc020/results/snpgenie/"
    conda:
        "../envs/env.yaml"
    output:
       "../results/snpgenie/{id}/population_summary.txt"
    shell:
        '''
        cd {params.dir}
        perl ../../workflow/scripts/snpgenie.pl \
            --minfreq=0.03 \
            --snpreport={input.vcf} \
            --vcfformat=4 \
            --fastafile={input.ref} \
            --gtffile={input.gtf} \
            --slidingwindow=150 \
            --outdir={wildcards.id}
        '''