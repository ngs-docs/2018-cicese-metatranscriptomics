import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version

from snakemake.remote import FTP
FTP = FTP.RemoteProvider()

min_version("5.1.2") #minimum snakemake version

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
data_targs, trim_targs, sourmash_targs = [], [], []
extensions = ["_1.fq.gz", "_2.fq.gz"]
trim_ext = ["_1.trim.fq.gz", "_2.trim.fq.gz", "_1.se.trim.fq.gz", "_2.se.trim.fq.gz"]
sourmash_ext =  ["_1.trim.sig", "_2.trim.sig"]

for s in SAMPLES: 
    data_targs = data_targs +  [s + e for e in extensions]
    trim_targs = trim_targs + [s + e for e in trim_ext]
    sourmash_targs = sourmash_targs + [s + e for e in sourmash_ext]

TARGETS = []

#def generate_data_targs(outdir, samples, extensions):
#    target_list = []
#    for s in samples:
#        target_list = target_list + [s + e for e in extensions]
        #target_list = target_list + [join(outdir, s + e) for e in extensions]
    #return target_list
#    return [join(outdir, t) for t in target_list]

## working here
#def generate assembly_targs(outdir, base, extensions):
#    target_list = target_list + [join(outdir, s + e) for e in extensions]

if download_data:
    include: join(RULES_DIR, 'general', 'ftp.rule')
else:
    include: join(RULES_DIR, 'general', 'link_data.rule')

include: join(RULES_DIR, 'trimmomatic', 'trimmomatic.rule')
#data_targs = [join(DATA_DIR, targ) for targ in data_targs] 
# Assembly Targets
include: join(RULES_DIR, 'trinity', 'trinity.rule')
trinity_extensions = ['_trinity.fasta', '_trinity.fasta.gene_trans_map']
trinity_targets = ['{}'.format(BASE) + i for i in trinity_extensions]
trinity_targs = [join(ASSEMBLY_DIR, t) for t in trinity_targets]

include: join(RULES_DIR, 'spades', 'spades.rule')
spades_targets = [join(ASSEMBLY_DIR, t) for t in [BASE + '_spades.fasta']]

include: join(RULES_DIR, 'plass', 'plass.rule')
plass_targets = [join(ASSEMBLY_DIR, t) for t in [BASE + '_plass.fasta']]

include: join(RULES_DIR, 'megahit', 'megahit.rule')
megahit_targets = [join(ASSEMBLY_DIR, t) for t in [BASE + '_megahit.fasta']]

#TARGETS = TARGETS + download_targs + [join(TRIM_DIR, targ) for targ in trim_targs] #+ trinity_targs
#TARGETS =  [join(TRIM_DIR, targ) for targ in trim_targs]
#TARGETS = spades_targets + plass_targets + megahit_targets #+ trinity_targs + sourmash_targets
#TARGETS = [join(TRIM_DIR, targ) for targ in trim_targs]

include: join(RULES_DIR, 'sourmash', 'sourmash.rule')
TARGETS = [join(TRIM_DIR, targ) for targ in sourmash_targs]

rule all:
    input: TARGETS


