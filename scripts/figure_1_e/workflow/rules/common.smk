import pandas as pd
import glob

def get_samples(wildcards):
    SAMPLES=pd.read_csv(config["sample_file"])['Sample']
    return SAMPLES

def get_r1(wildcards):
    #given ids, expand with glob to retrieve the full file name with added nucleotides
    return glob.glob(config["in_dir"]+"/"+wildcards.sample+'_*R1*.fastq.gz')

def get_r2(wildcards):
    return glob.glob(config["in_dir"]+"/"+wildcards.sample+'_*R2*.fastq.gz')
