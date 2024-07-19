rule index_vcf:
    input:
        A = "../results/snvs/{id}_A.vcf.gz",
        B = "../results/snvs/{id}_B.vcf.gz"
    output:
        A = "../results/snvs/{id}_A.vcf.gz.csi",
        B = "../results/snvs/{id}_B.vcf.gz.csi"
    shell:
        '''
        module load bcftools
        bcftools index -f {input.A} > {output.A}
        bcftools index -f {input.B} > {output.B}
        '''

rule filter_snvs:
    input:
        A_vcf = "../results/snvs/{id}_A.vcf.gz",
        B_vcf = "../results/snvs/{id}_B.vcf.gz",
        A_index = "../results/snvs/{id}_A.vcf.gz.csi",
        B_index = "../results/snvs/{id}_B.vcf.gz.csi"
    params:
        dir = "../results/filtered_snvs/{id}"
    output:
        expand("../results/filtered_snvs/{{id}}/{num}.vcf.gz", num = ["0000","0001","0002","0003"])
    shell:
        '''
        module load bcftools
        bcftools isec -c none -Oz {input.A_vcf} {input.B_vcf} -p {params.dir}
        touch {output}
        '''

rule merge_vcfs:
    input:
        expand("../results/filtered_snvs/{{id}}/{num}.vcf.gz", num = ["0002", "0003"])
    output:
        "../results/filtered_snvs/{id}/{id}_filtered_merged.vcf"
    shell:
        '''
        module load bcftools
        bcftools merge --merge all {input[0]} {input[1]} | bcftools +fill-tags > {output}
        '''