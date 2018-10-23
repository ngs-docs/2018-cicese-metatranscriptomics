import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version

from snakemake.remote import FTP
FTP = FTP.RemoteProvider()

configfile: "tara_config.yaml"

min_version("5.1.2") #minimum snakemake version

# read in sample info 
units = pd.read_table(config["units"]).set_index("sample", drop=False)

DATADIR = "/tara/data/TS_RNA"

samples = units['sample'].tolist()
targs = []
extensions = ["_1.fq.gz", "_2.fq.gz"]
for s in samples:
    targs = targs +  [s + e for e in extensions]

TARGETS = [join(DATADIR, targ) for targ in targs]

rule all:
    input: TARGETS

rule get_fq1:
    input: lambda wildcards: FTP.remote(expand("{file}", file=units.loc[wildcards.sample, "fq1"]), static=True, keep_local=True)
	output: join(DATADIR,"{sample}_1.fq.gz")
	shell: "cp {input} {output}"

rule get_fq2:
    input: lambda wildcards: FTP.remote(expand("{file}", file=units.loc[wildcards.sample, "fq2"]), static=True, keep_local=True)
	output: join(DATADIR,"{sample}_2.fq.gz")
	shell: "cp {input} {output}"

# use fastqc & trimmomatic rules from eelpond
# assemblers! 
  # plass rule
  # trinity rule -- use eelpond rule
  # megahit rule
  # rnaSpades rule
	
