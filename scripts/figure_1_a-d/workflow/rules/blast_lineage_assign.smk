rule concatenate_blast:
    input:
        expand("../results/blast_results/{exp}.{num}.out", exp=EXP, num=NUM)
    output:
        expand("../results/blast/{exp}_sorted_blast.txt", exp=EXP)
    shell:
        '''
        cat {input} | sort -k 1,1 -k4,4nr | sort -u -k1,1 --merge > {output}
        '''

rule assign_lineages:
    input:
        blast = expand("../results/blast/{exp}_sorted_blast.txt", exp=EXP),
    output:
       expand("../results/final/{exp}_transcript_lineages.csv", exp=EXP)
    params:
        nodes = config["blast"]["tax_nodes"],
        names = config["blast"]["tax_names"]
    shell:
        '''
        cut -f6,7 {input.blast} | sed 's/;.*//;s/\t/,/g' > ../results/blast/blast_lineage.csv
        python ~/bin/2018-ncbi-lineages/make-lineage-csv.py \
            {params.nodes} {params.names} \
            ../results/blast/blast_lineage.csv \
            -o {output}
        sed -i '1i trinity_id\tqseqid\tsseqid\tevalue\tbitscore\tsgi\tsacc\tstaxids\tsscinames\tscomnames\tstitle' {input.blast}
        cut -f1 {input.blast} | paste -d, - {output} > tmp.csv
        mv tmp.csv {output}
        rm ../results/blast/blast_lineage.csv
        '''