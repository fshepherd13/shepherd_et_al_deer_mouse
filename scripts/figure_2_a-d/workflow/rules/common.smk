import pandas as pd
import glob

def get_r1(wildcards):
    #given ids, expand with glob to retrieve the full file name with added nucleotides
    return glob.glob(config["in_dir"]+"/"+wildcards.id+"_"+wildcards.rep+wildcards.pool+'_*R1*.fastq.gz')

def get_r2(wildcards):
    return glob.glob(config["in_dir"]+"/"+wildcards.id+"_"+wildcards.rep+wildcards.pool+'_*R2*.fastq.gz')