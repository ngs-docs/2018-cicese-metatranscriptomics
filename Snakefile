import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version


configfile: "config.yaml"

min_version("5.1.2") #minimum snakemake version

# read in & validate sample info 
tara_samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

def get_links(wildcards):
    #return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    return dict(zip(['fq1','fq2'], units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()))

# download ftp

rule download_ftp:
    input:
	    unpack(get_links)
	output: 
	    fq1 = {sample}_1.fq.gz,
		fq2 = {sample}_2.fq.gz
    run:
        shell("wget {input.fq1} -O {output.fq1}")
        shell("wget {input.fq2} -O {output.fq2}")


get_ftp_targets(sampleInf):
	extensions = ['1.fq.gz', '2.fq.gz']
    samples = [sampleInf["sample"].tolist())
	targs = []
	for sample in samples:
        targs = targs +  ['{}_'.format(sample) + e for e in extensions]

TARGETS = get_ftp_targets(tara_samples)

# use fastqc & trimmomatic rules from eelpond
# assemblers! 
  # plass rule
  # trinity rule -- use eelpond rule
  # megahit rule
  # rnaSpades rule

rule all:
    input:
	    TARGETS

