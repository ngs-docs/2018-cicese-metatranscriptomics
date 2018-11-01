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
def is_single_end(sample, unit, end = ''):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def generate_data_targs(outdir, samples, extensions, ends = ["_1", "_2"]):
    target_list = []
    # to do: add paired vs single end check here to generate `ends`
    exts = [x+y for x in ends for y in extensions]
    for s in samples:
        target_list = target_list + [join(outdir, s + e) for e in exts]
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
QC_DIR = join(OUT_DIR, "read_qc")
ASSEMBLY_DIR = join(OUT_DIR,"assembly")
PALADIN_DIR = join(OUT_DIR,"paladin")

SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()

TARGETS = []

# download or softlink data
if download_data:
    include: join(RULES_DIR, 'general', 'ftp.rule')
else:
    include: join(RULES_DIR, 'general', 'link_data.rule')

data_ext = [".fq.gz", ".fq.gz"]
data_targs = generate_data_targs(DATA_DIR, SAMPLES, data_ext)

#trimmomatic trimming
include: join(RULES_DIR, 'trimmomatic', 'trimmomatic.rule')
trim_ext = [".trim.fq.gz", ".se.trim.fq.gz"]
trim_targs = generate_data_targs(TRIM_DIR, SAMPLES, trim_ext)

#fastqc of raw, trimmed files
#include: join(RULES_DIR,'fastqc/fastqc.rule')
#fastqc_pretrim_ext =  ['_fastqc.zip','_fastqc.html'] 
#fastqc_trim_ext =  ['_trimmed_fastqc.zip','_trimmed_fastqc.html'] 
#fastqc_pretrim_targs = [] #generate_data_targs(QC_DIR, SAMPLES, fastqc_pretrim_ext)
#fastqc_trim_targs = generate_data_targs(QC_DIR, SAMPLES, fastqc_trim_ext)
#fastqc_targs = fastqc_pretrim_targs + fastqc_trim_targs

#correct reads with rcorrector
include: join(RULES_DIR,'rcorrector/rcorrector.rule')
rcorr_ext =  ['.rcorr.fq.gz', '.se.rcorr.fq.gz']
rcorr_targs = generate_data_targs(TRIM_DIR, SAMPLES, rcorr_ext)

# error trim with khmer
include: join(RULES_DIR, 'khmer/khmer.rule')
khmer_pe_ext = ['_1.khmer.fq.gz', '_2.khmer.fq.gz', '.paired.khmer.fq.gz', '.single.khmer.fq.gz']
khmer_targs = generate_data_targs(TRIM_DIR, SAMPLES, khmer_pe_ext, ends = [""])

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

# generate sourmash signatures of trimmed reads and assemblies
include: join(RULES_DIR, 'sourmash', 'sourmash.rule')
sourmash_read_ext =  [".trim.sig"] 
sourmash_targs = generate_data_targs(TRIM_DIR, SAMPLES, sourmash_read_ext)
sourmash_assemb_ext = ['_megahit.sig', '_trinity.sig', '_plass.sig', '_spades.sig']
sourmash_targs = sourmash_targs + generate_base_targs(ASSEMBLY_DIR, BASE, sourmash_assemb_ext)

#include: join(RULES_DIR, 'paladin/paladin.rule')
#paladin_ext =  [".bwt"] 
#paladin_targs =  generate_base_targs(PALADIN_DIR, BASE, paladin_ext)

#TARGETS = TARGETS + download_targs + [join(TRIM_DIR, targ) for targ in trim_targs] #+ trinity_targs
#TARGETS =  [join(TRIM_DIR, targ) for targ in trim_targs]
#TARGETS = spades_targets + plass_targets + megahit_targets #+ trinity_targs + sourmash_targets
#TARGETS = [join(TRIM_DIR, targ) for targ in trim_targs]
#TARGETS = fastqc_targs + sourmash_targs
#TARGETS  =  khmer_targs
TARGETS =  trinity_targs # + spades_targs
#TARGETS = fastqc_targs

rule all:
    input: TARGETS


