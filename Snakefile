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
    target_list = [join(outdir, basename + e) for e in extensions]
    return target_list

# set up dirs, basename
BASE = config.get('basename','tara')
experiment_suffix = config.get('experiment_suffix')

if experiment_suffix:
    OUT_DIR = BASE + "_out_" + experiment_suffix
else:
    OUT_DIR = BASE + '_out'

RULES_DIR = 'rules'
DATA_DIR = config.get('data_directory', join(OUT_DIR, 'data'))
LOGS_DIR = join(OUT_DIR, 'logs')
TRIM_DIR = join(OUT_DIR,"trimmed")
QC_DIR = join(OUT_DIR, "read_qc")
ASSEMBLY_DIR = join(OUT_DIR,"assembly")
PALADIN_DIR = join(OUT_DIR,"paladin")
QUANT_DIR = join(OUT_DIR, 'quant')

# read in sample info 
samples = pd.read_table(config["samples"],dtype=str).set_index(["sample", "unit"], drop=False)
SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()

# determine workflow
flow = config.get('workflow', 'full')
read_processing,assembly,mapping,assembly_quality,annotation,quantification,diffexp,input_assembly = [False]*8

if flow == 'full':
    read_processing = True
    assembly = True
    assembly_quality = True
    annotation = True
    quantification = True
    diffexp = True
elif flow =='assembly':
#    read_processing = True
    mapping = True
    assembly = True
    quality = True
else:
    input_assembly = True
    assert config['assembly_input']['assembly'] is not None, "chosen workflow requires transcriptome assembly as input"
    assert config['assembly_input']['assembly_ext'] is not None, "chosen workflow requires transcriptome assembly as input"
    #assert config['assembly_input']['gene_trans_map'] is not None, "chosen workflow requires transcriptome gene transcript map as input"
    if flow == 'mapping':
        read_processing = True
        mapping = True
    #if flow == 'annotation':
    #    annotation = True
    #if flow == 'expression':
    #    read_processing = True
    #    quantification = True
    #    diffexp = True

TARGETS = []
include: 'rules/common.rule'

if read_processing:
    # download or softlink data
    if config.get('download_data', False):
        include: join(RULES_DIR, 'general', 'ftp.rule')
    else:
        include: join(RULES_DIR, 'general', 'link_data.rule')

    data_ext = [".fq.gz", ".fq.gz"]
    data_targs = generate_data_targs(DATA_DIR, SAMPLES, data_ext)
    #fastqc of raw, trimmed files
    #include: join(RULES_DIR,'fastqc/fastqc.rule')
    #fastqc_pretrim_ext =  ['_fastqc.zip','_fastqc.html'] 
    #fastqc_trim_ext =  ['_trimmed_fastqc.zip','_trimmed_fastqc.html'] 
    #fastqc_pretrim_targs = [] #generate_data_targs(QC_DIR, SAMPLES, fastqc_pretrim_ext)
    #fastqc_trim_targs = generate_data_targs(QC_DIR, SAMPLES, fastqc_trim_ext)
    #fastqc_targs = fastqc_pretrim_targs + fastqc_trim_targs

    #trimmomatic trimming
    include: join(RULES_DIR, 'trimmomatic', 'trimmomatic.rule')
    trim_ext = [".trim.fq.gz", ".se.trim.fq.gz"]
    trim_targs = generate_data_targs(TRIM_DIR, SAMPLES, trim_ext)

    #correct reads with rcorrector
    include: join(RULES_DIR,'rcorrector/rcorrector.rule')
    rcorr_ext =  ['.rcorr.fq.gz']
    rcorr_targs = generate_data_targs(TRIM_DIR, SAMPLES, rcorr_ext)

    # error trim with khmer
    include: join(RULES_DIR, 'khmer/khmer.rule')
    khmer_pe_ext = ['_1.khmer.fq.gz', '_2.khmer.fq.gz', '.paired.khmer.fq.gz', '.single.khmer.fq.gz']
    khmer_targs = generate_data_targs(TRIM_DIR, SAMPLES, khmer_pe_ext, ends = [""])

    # generate sourmash signatures of trimmed reads
    include: join(RULES_DIR, 'sourmash', 'sourmash.rule')
    sourmash_read_ext =  [".trim.sig"] 
    sourmash_targs = generate_data_targs(TRIM_DIR, SAMPLES, sourmash_read_ext)

if assembly:
    assemblies = []
    # trinity assembly
    include: join(RULES_DIR, 'trinity', 'trinity.rule')
    trinity_ext = ['_trinity.fasta', '_trinity.fasta.gene_trans_map']
    trinity_targs = generate_base_targs(ASSEMBLY_DIR, BASE, trinity_ext)
    assemblies+=['trinity']
    
    # spades assembly
    include: join(RULES_DIR, 'spades', 'spades.rule')
    spades_ext = ['_spades.fasta']
    spades_targs = generate_base_targs(ASSEMBLY_DIR, BASE, spades_ext)
    assemblies+=['spades']

    # plass assembly
    include: join(RULES_DIR, 'plass', 'plass.rule')
    plass_ext = ['_plass.fasta']
    plass_targs = generate_base_targs(ASSEMBLY_DIR, BASE, plass_ext)
    assemblies+=['plass']

    # megahit assembly
    include: join(RULES_DIR, 'megahit', 'megahit.rule')
    megahit_ext = ['_megahit.fasta']
    megahit_targs = generate_base_targs(ASSEMBLY_DIR, BASE, megahit_ext)
    assemblies+=['megahit']

    # generate sourmash signatures of assemblies
    include: join(RULES_DIR, 'sourmash', 'sourmash.rule')
    #sourmash_assemb_ext = ['_megahit.sig', '_trinity.sig', '_plass.sig', '_spades.sig']
    sourmash_assemb_ext = ['_' + x + '.sig' for x in assemblies]
    sourmash_assemb_targs = generate_base_targs(ASSEMBLY_DIR, BASE, sourmash_assemb_ext)

if input_assembly:
    include: 'rules/assemblyinput/assemblyinput.rule'
    assemb = BASE + str(config['assembly_input']['assembly_ext'])
    assemblyinput_targs = generate_base_targs(ASSEMBLY_DIR, assemb, ['.fasta'])

if mapping:
    #PEAR: merge pe reads for paladin
#    include: join(RULES_DIR, 'pear/pear.rule')
#    pear_ext = ['.pear_merged.fq.gz']
#    pear_targs = generate_data_targs(TRIM_DIR, SAMPLES, pear_ext, ends = [""])
    #paladin_read_ext =  ["_pear.paladin.bam", "_pear.paladin.sort.bam", "_pear.paladin.sort.bam.bai"] 
    # PALADIN: map reads in AA space
    include: join(RULES_DIR, 'paladin/paladin.rule')
    paladin_read_ext =  ["_trim.paladin.bam", "_trim.paladin.sort.bam", "_trim.paladin.sort.bam.flagstat" , "_trim.paladin.sort.bam.bai"] 
    assemb_name =BASE + '_plass'
    paladin_targs =  generate_base_targs(PALADIN_DIR + '_' + assemb_name, assemb_name, [".fasta", ".fasta.bwt"])
    paladin_targs += generate_data_targs(PALADIN_DIR + '_' + assemb_name , SAMPLES, paladin_read_ext)
    include: join(RULES_DIR, 'salmon/salmon.rule')
    salmon_read_ext = ['/quant.sf', '/lib_format_counts.json']
    salmon_targs = []
    for assembly in assemblies:
        if assembly != 'plass':
            salmon_targs += generate_base_targs(join(QUANT_DIR, BASE + '_' + assembly), BASE + '_' + assembly, ['_salmon']) 
            salmon_targs += generate_data_targs(join(QUANT_DIR, BASE + '_' + assembly), SAMPLES, salmon_read_ext, ends = [""])

#TARGETS = TARGETS + download_targs + [join(TRIM_DIR, targ) for targ in trim_targs] #+ trinity_targs
#TARGETS =  [join(TRIM_DIR, targ) for targ in trim_targs]
#TARGETS = spades_targets + plass_targets + megahit_targets #+ trinity_targs + sourmash_targets
#TARGETS = [join(TRIM_DIR, targ) for targ in trim_targs]
#TARGETS = fastqc_targs + sourmash_targs
#TARGETS  =  khmer_targs
#TARGETS =  trinity_targs # + spades_targs
#TARGETS = fastqc_targs
#TARGETS = paladin_targs + salmon_read_targs
#TARGETS = rcorr_targs 
#TARGETS = megahit_targs
#TARGETS = plass_targs + megahit_targs
TARGETS = trinity_targs + salmon_targs + paladin_targs

rule all:
    input: TARGETS


