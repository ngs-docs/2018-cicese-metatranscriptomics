"""Snakemake wrapper for rnaSpades."""

__author__ = "Tessa Pierce"
__copyright__ = "Copyright 2018, Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
memory =  snakemake.params.get("memory", "64G")

#allow multiple input files for single assembly
left = snakemake.input.get("left")
if left:
    right = snakemake.input.get("right")

single = snakemake.input.get("single")
interleaved = snakemake.input.get("interleaved")
assert left is not None or single is not None or interleaved is not None, "please check read input"

input_cmd = ''

if left:
    left = [snakemake.input.left] if isinstance(snakemake.input.left, str) else snakemake.input.left
    right = [snakemake.input.right] if isinstance(snakemake.input.right, str) else snakemake.input.right 
    assert len(left) == len(right), "left input needs to contain the same number of files as the right input" 
    input_str_left = ' -1 ' + " -1 ".join(left)
    input_str_right = ' -2 ' + " -2 ".join(right)
    input_cmd =  input_cmd + " ".join([input_str_left, input_str_right])
if interleaved:
    interleaved = [snakemake.input.interleaved] if isinstance(snakemake.input.interleaved, str) else snakemake.input.interleaved
    input_str_interleaved = ' --12 ' + " --12 ".join(interleaved)
    input_cmd =  input_cmd + " ".join([input_str_interleaved])
if single:
    single = [snakemake.input.single] if isinstance(snakemake.input.single, str) else snakemake.input.single 
    input_str_single = ' -s ' + " -s ".join(single)
    input_cmd =  input_cmd + " ".join([input_str_single])

all_files = left + right + interleaved + single
assert len(all_files) <=9, "rnaspades can only take 9 input files via the command line interface"

outdir = path.dirname(snakemake.output[0])

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

#shell("Trinity {input_cmd} --CPU {snakemake.threads} --max_memory {max_memory} --seqType {seqtype} --output {outdir} {snakemake.params.extra} {log}")
shell("spades.py --rna -o {outdir} {input_cmd} -t {snakemake.threads} -m {memory} {extra}") 

