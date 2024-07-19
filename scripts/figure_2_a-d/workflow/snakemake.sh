#Part 1: Get consensus
snakemake -s Snakefile_get_consensus --cores 54 --rerun-incomplete --use-conda --configfile ../config/config_aj34.yaml

#Part 2: Check primers
snakemake --cores 54 --rerun-incomplete -s Snakefile_check_primers --use-conda --configfile ../config/config_aj34.yaml

#Part 2: Call variants
snakemake --cores 54 --rerun-incomplete -s Snakefile_get_variants --use-conda --configfile ../config/config_aj34.yaml
