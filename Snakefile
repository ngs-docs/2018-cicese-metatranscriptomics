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
samples = pd.read_table(config["samples"]).set_index(["sample", "unit"], drop=False)

BASE = config.get('basename','tara')
experiment_suffix = config.get('experiment_suffix', '')

OUT_DIR = BASE + "_" + experiment_suffix

LOGS_DIR = join(OUT_DIR, 'logs')
DATA_DIR = join(OUT_DIR, "/tara/data/TS_RNA")
TRIM_DIR = join(OUT_DIR,"trimmed")
ASSEMBLY_DIR = join(OUT_DIR,"assembly")

#get ftp download targets
SAMPLES = (samples['sample'].astype(str) + '_' + samples['unit'].astype(str)).tolist()
ftp_targs, trim_targs = [], []
extensions = ["_1.fq.gz", "_2.fq.gz"]
trim_ext = ["_1.trim.fq.gz", "_2.trim.fq.gz", "_1.se.trim.fq.gz", "_2.se.trim.fq.gz"]
for s in SAMPLES: 
    ftp_targs = ftp_targs +  [s + e for e in extensions]
    trim_targs = trim_targs + [s + e for e in trim_ext]

TARGETS = []

# if some of the files are missing, add ftp targets (aka download the files!)
ftp_targs = [join(DATA_DIR, targ) for targ in ftp_targs] 
download_targs = []
if not all([os.path.isfile(f) for f in ftp_targs]):  
    for f in ftp_targs:
        if not os.path.isfile(f):
            download_targs = download_targs + [f] # only download missing files  

# Assembly Targets
trinity_extensions = ['_trinity.fasta', '_trinity.fasta.gene_trans_map']
trinity_targets = ['{}'.format(BASE) + i for i in trinity_extensions]
trinity_targs = [join(ASSEMBLY_DIR, t) for t in trinity_targets]


#TARGETS = TARGETS + download_targs + [join(TRIM_DIR, targ) for targ in trim_targs] #+ trinity_targs
#TARGETS =  [join(TRIM_DIR, targ) for targ in trim_targs]
#TARGETS = trinity_targs
TARGETS = [join(TRIM_DIR, targ) for targ in trim_targs]


#find files for fastqc & trimmomatic input
def get_pretrim_pe(wildcards):
    # if specifying file locations in csv, rather than download links, use this
    #return dict(zip(['r1','r2'], samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()))
    # if download links are specified, we need to grab the downloaded files instead
    fq1 = join(DATA_DIR, '{}_{}_1.fq.gz'.format(wildcards.sample,wildcards.unit))
    fq2 = join(DATA_DIR, '{}_{}_2.fq.gz'.format(wildcards.sample,wildcards.unit))
    return dict(zip(['r1','r2'],[fq1,fq2]))

rule all:
    input: TARGETS

rule get_fq1:
    input: lambda wildcards: FTP.remote(expand("{file}", file=samples.loc[[wildcards.sample,wildcards.unit], "fq1"]), static=True, keep_local=True, immediate_close=True)
    output: join(DATA_DIR,"{sample}_{unit}_1.fq.gz")
    shell: "mv {input} {output}"

rule get_fq2:
    input: lambda wildcards: FTP.remote(expand("{file}", file=samples.loc[[wildcards.sample,wildcards.unit], "fq2"]), static=True, keep_local=True, immediate_close=True)
    output: join(DATA_DIR,"{sample}_{unit}_2.fq.gz")
    shell: "mv {input} {output}"


trim_params = config['trimmomatic']

rule trimmomatic_pe:
    """
    Trim reads from the sequencer by trimming or dropping low-quality reads.
    """
    input:
	    r1= lambda wildcards: join(DATA_DIR, '{}_{}_1.fq.gz'.format(wildcards.sample,wildcards.unit)), #unpack(get_pretrim_pe) 
	    r2= lambda wildcards: join(DATA_DIR, '{}_{}_2.fq.gz'.format(wildcards.sample,wildcards.unit))

    output:
        r1=join(TRIM_DIR, "{sample}_{unit}_1.trim.fq.gz"),
        r2=join(TRIM_DIR, "{sample}_{unit}_2.trim.fq.gz"),
        r1_unpaired=join(TRIM_DIR, "{sample}_{unit}_1.se.trim.fq.gz"),
        r2_unpaired=join(TRIM_DIR, "{sample}_{unit}_2.se.trim.fq.gz"),
    message:
        """--- Quality trimming PE read data with Trimmomatic."""
    params:
        trimmer = (trim_params['trim_cmd'].format(trim_params['adapter_file']['pe_name'])).split(' '),
        extra = '' 
    log: join(LOGS_DIR, 'trimmomatic/{sample}_{unit}_pe.log')
    conda:"trimmomatic-env.yaml"
    script:"trimmomatic-pe.py"

rule trinity:
    input:
            left=expand(join(TRIM_DIR, '{sample}_{end}.trim.fq.gz'), sample=SAMPLES, end=["1", "1.se","2.se"]), 
	    right=expand(join(TRIM_DIR, '{sample}_2.trim.fq.gz'), sample=SAMPLES)
    output:
        fasta = join(ASSEMBLY_DIR,"trinity_out_dir/Trinity.fasta"),
        gene_trans_map = join(ASSEMBLY_DIR,"trinity_out_dir/Trinity.fasta.gene_trans_map"),
    message:
        """### Assembling read data with Trinity ### """
    params:
        #**config['trinity']
        # optional parameters
        seqtype='fq',
        max_memory='64G',
        extra=""
    threads: 4
    log: join(LOGS_DIR, 'trinity/trinity.log')
    conda: "trinity-env.yaml"
    script: "trinity-wrapper.py"

rule rename_trinity_fasta:
    input: rules.trinity.output.fasta
    output: join(ASSEMBLY_DIR, BASE + '_trinity.fasta')
    log: join(LOGS_DIR, 'trinity/cp_assembly.log')
    shell: ("cp {input} {output}") 

rule rename_trinity_gene_trans_map:
    input: rules.trinity.output.gene_trans_map
    output: join(ASSEMBLY_DIR, BASE + '_trinity.fasta.gene_trans_map')
    log: join(LOGS_DIR, 'trinity/cp_gt_map.log')
    shell: ("cp {input} {output}") 


#rule rnaspades:
#    input:
#            left=expand(join(TRIM_DIR, '{sample}_{end}.trim.fq.gz'), sample=SAMPLES, end=["1", "1.se","2.se"]), 
#	    right=expand(join(TRIM_DIR, '{sample}_2.trim.fq.gz'), sample=SAMPLES)
#    output:
#        fasta = join(ASSEMBLY_DIR,"trinity_out_dir/Trinity.fasta"),
#        gene_trans_map = join(ASSEMBLY_DIR,"trinity_out_dir/Trinity.fasta.gene_trans_map"),
#    message:
#        """### Assembling read data with Trinity ### """
#    params:




# use fastqc & trimmomatic rules from eelpond
# assemblers! 
  # plass rule
  # megahit rule
  # rnaSpades rule
#	
