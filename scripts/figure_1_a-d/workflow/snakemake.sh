#AJ30 analysis:
snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj30/config_dm_aj30.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj30/config_mus_aj30.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part2 --configfile ../config/aj30/config_aj30.yaml

#AJ32 analysis:
snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj32/config_dm_aj30.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj32/config_mus_aj30.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part2 --configfile ../config/aj32/config_aj30.yaml

#AJ34 analysis:
snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj34/config_dm_aj34.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj34/config_mus_aj34.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part2 --configfile ../config/aj34/config_aj34.yaml

#AJ35 analysis:
snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj35/config_dm_aj35.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part1 --configfile ../config/aj35/config_mus_aj35.yaml

snakemake --cores 50 --rerun-incomplete --use-conda -s Snakefile_part2 --configfile ../config/aj35/config_aj35.yaml
