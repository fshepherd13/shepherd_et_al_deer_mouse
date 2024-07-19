#AJ39 analysis:
snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj39/config_dm_aj39.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj39/config_mus_aj39.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part2 --configfile ../config/aj39/config_aj39.yaml
