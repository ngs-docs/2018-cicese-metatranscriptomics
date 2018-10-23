import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version

#configfile: "tara_config.yaml"

min_version("5.1.2") #minimum snakemake version

# read in & validate sample info 
#units = pd.read_table(config["units"]).set_index("sample", drop=False)
units = pd.read_table('taraTSRNA_units.tsv', dtype=str).set_index(["sample"], drop=False)


# download ftp
#DATADIR = "/tara/data/TS_RNA"
DATADIR = ""

#def get_ftp_targets(sampleInf):
#    extensions = ["_1.fq.gz", "_2.fq.gz"]
#    samples = sampleInf["sample"].tolist()
    #print(samples)
#    targs = []
#    for s in samples:
#        print(s)
     #   fq1 = s + "_1.fq.gz"
     #   fq2 = s + "_2.fq.gz"
     #   import pdb;pdb.set_trace()
     #   targs = targs + [fq1,fq2]
#        targs = targs +  [s + e for e in extensions]
#	return targs

samples = units['sample'].tolist()
targs = []
extensions = ["_1.fq.gz", "_2.fq.gz"]
for s in samples:
    targs = targs +  [s + e for e in extensions]
#targs= get_ftp_targets(units)


TARGETS = [join(DATADIR, targ) for targ in targs]

rule all:
    input: TARGETS

def get_links(wildcards):
    #return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    return dict(zip(["fq1","fq2"], units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()))

rule download_ftp:
    input:
        unpack(get_links)
    output: 
        fq1=join(DATADIR,"{sample}_1.fq.gz"),
        fq2=join(DATADIR,"{sample}_2.fq.gz")
    run:
        shell("wget {input.fq1} -O {output.fq1}")
        shell("wget {input.fq2} -O {output.fq2}")

rule trimmomatic:
    input:
       


# use fastqc & trimmomatic rules from eelpond
# assemblers! 
  # plass rule
  # trinity rule -- use eelpond rule
  # megahit rule
  # rnaSpades rule
	
