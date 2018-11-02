__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2018, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
k_size = snakemake.params.get("k", "31")

r1 = snakemake.input.get("r1")
r2 = snakemake.input.get("r2")
#r1_out = snakemake.output.get("r1")
#r2_out = snakemake.output.get("r2")

ileav = snakemake.input.get("i")
#ileav_out = snakemake.output.get("i")

unp = snakemake.input.get('s')
#unp_out = snakemake.output.get("s")


# maybe generate output tuples here: (default_outfile, desired_outfile)
def generate_default_output(infiles, outD):
    outfiles = []
    for f_in in infiles:
        base = path.basename(f_in.rsplit('.f')[0])
        default_outname = base + '.rcor.fq'
        # Rcorrector outputs gzipped files IF input files are gzipped
        if f_in.endswith('.gz'):
            default_outname = default_outname + '.gz'
        outfiles = outfiles + [default_outname]
    return [path.join(outD, f) for f in outfiles]

assert (r1 is not None and r2 is not None) or ileav is not None or unp is not None, "either r1 and r2 (paired), or i (interleaved) are required as input"

if r1:
    # make sure we're handling list version of input/output
    r1 = [snakemake.input.r1] if isinstance(snakemake.input.r1, str) else snakemake.input.r1
    r2 = [snakemake.input.r2] if isinstance(snakemake.input.r2, str) else snakemake.input.r2
    r1_out = [snakemake.output.r1] if isinstance(snakemake.output.r1, str) else snakemake.output.r1
    r2_out = [snakemake.output.r2] if isinstance(snakemake.output.r2, str) else snakemake.output.r2
    #set outdir
    outdir = path.dirname(r1[0])
    # check number of input/output files
    assert len(r1) == len(r2), "input-> equal number of files required for r1 and r2"
    assert len(r1) == len(r1_out), "requires equal # output files as input files"
    assert len(r1) == len(r2_out), "requires equal # output files as input files"
    # generate command line for read input
    r1_cmd = ' -1 ' + " -1 ".join(r1)
    r2_cmd = ' -2 ' + " -2 ".join(r2)
    read_cmd = " ".join([r1_cmd,r2_cmd])
    # generate intermediate output names
    r1_default_out = generate_default_output(r1, outdir)
    r2_default_out = generate_default_output(r2, outdir)
if ileav:
    assert r1 is None and r2 is None and unp is None, "cannot mix paired, unpaired, and interleaved input" #to do: check this!
    ileav = [snakemake.input.i] if isinstance(snakemake.input.i, str) else snakemake.input.i
    ileav_out = [snakemake.output.i] if isinstance(snakemake.output.i, str) else snakemake.output.i
    #set outdir
    outdir = path.dirname(ileav[0])
    read_cmd = ' -i ' + " -i ".join(ileav)
    assert len(r1) == len(ileav_out), "requires equal # output files as input files"
    ileav_default_out = generate_default_output(ileav, outdir)
if unp:
    assert r1 is None and r2 is None and ileav is None, "cannot mix paired, unpaired, and interleaved input" #to do: check this!
    unp = [snakemake.input.s] if isinstance(snakemake.input.s, str) else snakemake.input.s
    unp_out = [snakemake.output.s] if isinstance(snakemake.output.s, str) else snakemake.output.s
    #set outdir
    outdir = path.dirname(unp[0])
    read_cmd = ' -s ' + " -s ".join(unp)
    assert len(unp) == len(unp_out), "requires equal # output files as input files"
    unp_default_out = generate_default_output(unp, outdir)

# run rcorrector
shell("run_rcorrector.pl {read_cmd} -od {outdir} -k {k_size} -t {snakemake.threads} {snakemake.params.extra} {log}")

# move default output to desired output name and location
if r1_default_out:
    for f_default, f_out in zip(r1_default_out, r1_out):
        shell("cp {f_default} {f_out}")
        shell("rm {f_default}")
    for f_default, f_out in zip(r2_default_out, r2_out):
        shell("cp {f_default} {f_out}")
        shell("rm {f_default}")

if ileav_default_out:
    for f_default, f_out in zip(ileav_default_out, ileav_out):
        shell("cp {f_default} {f_out}")
        shell("rm {f_default}")

if unp_default_out:
    for f_default, f_out in zip(unp_default_out, unp_out):
        shell("cp {f_default} {f_out}")
        shell("rm {f_default}")
