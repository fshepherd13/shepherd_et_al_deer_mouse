# Shepherd et al. Deer mouse manuscript code

This repository contains code to reproduce the analyses found in the deer mouse manuscript.

## Figure 1

### A-D

To rerun the RNAseq analysis, use the code found in `scripts/figure_1_a-d/`. Several exposure experiments were completed and therefore analyzed separately to generate estimated pathogen transcript counts. See Supplemental table S1 for experimental layouts and a list of each RNAseq analysis group. 

#### Dependencies:
To run the snakemake workflows, you will need:
1. An indexed reference genome for both Mus musculus and Peromyscus maniculatus
2. A blast database 
3. Taxonomy nodes and names file formatted for extracting taxonomic lineages from accession numbers. See the instructions at https://github.com/dib-lab/2018-ncbi-lineages
4. The code from https://github.com/dib-lab/2018-ncbi-lineages. Download and make sure the path in `rules/blast_lineage_assign.smk` rule `assign_lineages` references the correct path of the script. 

The Snakemake pipeline builds conda environments for each appropriate rule.

#### To run:

A snakemake pipeline to analyze the raw fastq files should be run in three parts. The first two parts are mapping the reads to the appropriate host genome. The last part is de novo assembling the remaining reads, identifying their taxonomic lineages, and quantifying pathogen reads per host. These steps can be achieved using the commands in `scripts/figure_1_a-d/workflow/snakemake.sh`. Each experiment (i.e. AJ30, AJ32, AJ34, and AJ35) is analyzed separately to most accurately de novo assemble the pathogen transcriptomes. This is because each experiment utilized a different pet store mouse, and their microbiomes/viromes can vary quite a bit. To avoid creating erroneous contigs that are combinations of different pathogen strains, only animals associated within an experimental batch were analyzed together. 

The outputs of the pipeline used for downstream analysis are Salmon estimated counts for each Trinity contig, and the lineages for each Trinity contig. For convenience, the output files are located in `scripts/figure_1_a-d/data`. These files are used as input for the Rscript found at `scripts/figure_1_a-d/Figure_1_a-d.Rmd` which normalizes read counts using DESeq2's median of ratios methods, and produces graphs for figure 1.

### E

The RNAseq analysis works the same for figure E as described above. Use the snakemake workflow within `scripts/figure_1_e/workflow` for reproducing the estimated pathogen counts with Salmon. Use the code in `scripts/figure_1_e/Figure_1_e.Rmd` to create the figure in panel E. 

## Figure 2
### A-D
Intra-host single nucleotide variants of Murine kobuvirus (MKV) were called using [iVar](https://github.com/andersen-lab/ivar) in a snakemake pipeline. Briefly, a consensus sequence was called using the day 0 pet store mouse feces reads mapped to a de novo assembled contig of MKV. The reads from all animals/samples are then mapped back to that consensus to call variants. 

To reproduce the analysis and variant calling, use the snakemake commands in order in the `scripts/figure_2_a-d/workflow/snakemake.sh` file. The steps were run separately in order to manually check for correct primer binding positions and for annotating the final consensus sequence output from the d0 feces reads. The final filtered tsv files are used as input for the R script, and they are included in this repo as they would appear after running the pipeline in `scripts/figure_2_a-d/results/variant_calling_AJ34_6_11_d0_FECES_MKV1_WG/final/`.

#### Dependencies
The Snakemake pipeline downloads the dependencies as conda environments for each rule. See the yaml files in `scripts/figure_2_a-d/workflow/envs`.

#### To run:

1. Step 1: get consensus. 
In this step, the reads from both replicates of the day 0 pet store mouse feces reads are merged and mapped to a reference MKV contig. The output is a fasta file of the consensus sequence which is used to manually annotated the MKV open reading frame for protein calling. Reads from all samples are also quality trimmed with bbduk and fastqc is run on all samples. The resulting consensus seq and gtf file is included in the repo.
2. Step 2: Check primers.
Here, the primers are mapped to the called consensus sequence to check for mismatches. Reads associated with amplicons generated by mismatched primers are discarded. 
3. Step 3: Call variants.
The trimmed reads generated in step 1 are mapped back to the annotated reference sequence and iSNVs are called with iVar using a minimum depth of 100x and minimum quality threshold of q30. 


### E

Figure E was created using [SNPGenie](https://github.com/chasewnelson/SNPGenie) and samtools/bcftools. The approach was as follows:

1. Get a BAM file of reads mapped to the d0 feces consensus sequence, one for each sequencing replicate (from the pipeline results from figure 1A-D above)
2. Generate pileup with samtools mpileup of the BAM files in VCF format and call variants with bcftools (rule call_snvs.smk)
3. Index the vcf files with bcftools index (rule filter_snvs.smk)
4. Filter variants with bcftools isec; this creates an output vcf file which lists the SNVs present in both sequendcing replicate BAM/pileup files (rule filter_snvs.smk).
5. Merge the vcf files so that there is one VCF file per animal containing the position, reference and alternate nucleotide information for the shared SNVs and a combined depth/frequency calculation for both replicates.
6. Edit the vcf headers with bcftools header (rule snpgenie.smk).
7. Perform piN/piS calculations with SNPgenie.

Steps 1-6 are done with a Snakemake pipeline, step 7 is done by hand (couldn't get snpgenie to work with the snakemake pipeline).

The snakemake pipeline uses a single environment which can be reproduced using the yaml file at `scripts/figure_2_e/config/env.yaml`. Then run the snakemake pipeline with `snakemake --configfile ../config/config_aj34.yaml`.

The output files are used as input to snpgenie, i.e (from within the `scripts/figure_2_e/results/snpgenie`):

```
perl snpgenie.pl \
    --minfreq=0.03 \
    --snpreport=AJ34_2_2_SI_MKV1_WG_filtered_merged_ed.vcf \
    --vcfformat=4 \
    --fastafile=../../../figure_2_a-d/results/consensus_calling/final/AJ34_6_11_d0_FECES_MKV1_WG.fa \
    --gtffile=../../../figure_2_a-d/results/consensus_calling/final/AJ34_6_11_d0_FECES_MKV1_WG.gtf \
    --slidingwindow=50 \
    --outdir=AJ34_2_2_SI_MKV1_WG
```

SNPgenie outputs multiple files but the `product_results.txt` file is what is used for publication/graphing. See the file at `scripts/figure_2_e/Figure_2_e.Rmd` for graphing code.