import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version

from snakemake.remote import FTP
FTP = FTP.RemoteProvider()

min_version("5.1.2") #minimum snakemake version

#util functions (to do: move elsewhere)
def generate_data_targs(outdir, samples, extensions):
    target_list = []
    for s in samples:
        target_list = target_list + [join(outdir, s + e) for e in extensions]
    return target_list

def generate_base_targs(outdir, basename, extensions):
    target_list = []
    target_list = [join(outdir, BASE + e) for e in extensions]
    return target_list

# read in sample info 
samples = pd.read_table(config["samples"],dtype=str).set_index(["sample", "unit"], drop=False)

BASE = config.get('basename','tara')
experiment_suffix = config.get('experiment_suffix')

if experiment_suffix:
    OUT_DIR = BASE + "_out_" + experiment_suffix
else:
    OUT_DIR = BASE + '_out'

RULES_DIR = 'rules'
DATA_DIR = config.get('data_directory', join(OUT_DIR, 'data'))
download_data = config.get('download_data', False)

LOGS_DIR = join(OUT_DIR, 'logs')
TRIM_DIR = join(OUT_DIR,"trimmed")
ASSEMBLY_DIR = join(OUT_DIR,"assembly")

SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()

TARGETS = []

# download or softlink data
if download_data:
    include: join(RULES_DIR, 'general', 'ftp.rule')
else:
    include: join(RULES_DIR, 'general', 'link_data.rule')

data_ext = ["_1.fq.gz", "_2.fq.gz"]
data_targs = generate_data_targs(DATA_DIR, SAMPLES, data_ext)

#trimmomatic trimming
include: join(RULES_DIR, 'trimmomatic', 'trimmomatic.rule')
trim_ext = ["_1.trim.fq.gz", "_2.trim.fq.gz", "_1.se.trim.fq.gz", "_2.se.trim.fq.gz"]
trim_targs = generate_data_targs(TRIM_DIR, SAMPLES, trim_ext)

# trinity assembly
include: join(RULES_DIR, 'trinity', 'trinity.rule')
trinity_ext = ['_trinity.fasta', '_trinity.fasta.gene_trans_map']
trinity_targs = generate_base_targs(ASSEMBLY_DIR, BASE, trinity_ext)

# spades assembly
include: join(RULES_DIR, 'spades', 'spades.rule')
spades_ext = ['_spades.fasta']
spades_targs = generate_base_targs(ASSEMBLY_DIR, BASE, spades_ext)

# plass assembly
include: join(RULES_DIR, 'plass', 'plass.rule')
plass_ext = ['_plass.fasta']
plass_targs = generate_base_targs(ASSEMBLY_DIR, BASE, plass_ext)

# megahit assembly
include: join(RULES_DIR, 'megahit', 'megahit.rule')
megahit_ext = ['_megahit.fasta']
megahit_targs = generate_base_targs(ASSEMBLY_DIR, BASE, megahit_ext)

#TARGETS = TARGETS + download_targs + [join(TRIM_DIR, targ) for targ in trim_targs] #+ trinity_targs
#TARGETS =  [join(TRIM_DIR, targ) for targ in trim_targs]
#TARGETS = spades_targets + plass_targets + megahit_targets #+ trinity_targs + sourmash_targets
#TARGETS = [join(TRIM_DIR, targ) for targ in trim_targs]

# generate sourmash signatures of trimmed reads and assemblies
include: join(RULES_DIR, 'sourmash', 'sourmash.rule')
sourmash_read_ext =  ["_1.trim.sig", "_2.trim.sig"] 
sourmash_targs = generate_data_targs(TRIM_DIR, SAMPLES, sourmash_read_ext)
sourmash_assemb_ext = ['_megahit.sig', '_trinity.sig', '_plass.sig', '_spades.sig']
sourmash_targs = sourmash_targs + generate_base_targs(ASSEMBLY_DIR, BASE, sourmash_assemb_ext)

TARGETS = sourmash_targs

rule all:
    input: TARGETS


