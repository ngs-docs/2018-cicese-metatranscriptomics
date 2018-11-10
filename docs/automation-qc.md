# Automation via shell scripts and snakemake

For everything we have done so far, we have copied and pasted a lot of commands 
to accomplish what we want. This works! But can also be time consuming, and is
more prone to error. We will show you next how to put all of these commands into
a shell script.
  
A **shell script** is a text file full of shell commands, that run just as if you're 
running them interactively at the command line.

Here, we'll create one that takes the commands from the and runs 'em all at once.

## Writing a shell script

Let's put all of our commands from the quality control lesson into one script. 

We'll call it `run_qc.sh`. The `sh` at the end of the tells you that this is a bash script. 

Using nano, edit the file `run-qc.sh` and put the following content there:

```
cd ${PROJECT}
mkdir -p quality
cd quality

ln -s ../data/*.fq.gz ./

printf "I see $(ls -1 *.fq.gz | wc -l) files here.\n"

fastqc *.fq.gz

multiqc .
```

This is now a shell script that you can use to execute all
of those commands in *one* go -- try it out! Run:

```
cd ~/
bash run-qc.sh
```

### Re-running the shell script

Suppose you wanted to re-run the script. How would you do that?

Well, note that the `quality` directory is created at the top of the script, and everything is executed in that directory. So if you remove the quality directory like so,

```
rm -fr quality
```

you can then do

```
bash run-qc.sh
```

### Some tricks for writing shell scripts

#### Make it executable

You can get rid of the `bash` part of the command above with
some magic:

Put 
```
#! /bin/bash
```
at the top of the file, and then run

```
chmod +x ~/run-qc.sh
```

at the command line.

You can now run

```
./run-qc.sh
```
instead of `bash run-rnaseq.sh`.

You might be thinking, ok, why is this important? Well, you can do the same with R scripts and Python scripts (but put `/usr/bin/env Rscript` or `/usr/bin/env python` at the top, instead of `/bin/bash`). This basically annotates the script with the language it's written in, so you don't have to know or remember yourself.

So: it's not necessary but it's a nice trick.

You can also always *force* a script to be run in a particular language by specifying `bash <scriptname>` or `Rscript <Scriptname>`, too.

#### Another good addition: make it fail nicely

One sort of weird aspect of shell scripts is that (by default) they keep running even if there's an error.  This is bad behavior and we should turn it off.

You can observe this behavior by rerunning the script above without deleting the directory `rnaseq/` - the `mkdir` command will print an error because the directory still exists, but 

A good addition to every shell script is to make it fail on the first error. Do this by putting
```
set -e
``` 
at the top - this tells bash to exit at the first error, rather than continuing bravely onwards.

#### A final good addition: make shell scripts print out the commands they're running!

You might notice that the shell script gives you the output of the commands it's running, but doesn't tell you *which* commands it's running.

If you add
```
set -x
```
at the top of the shell script and then re-run it,

```
cd ~/
rm -fr quality
./run-qc.sh
```

then you'll see the full set of commands being run!

#### A final note on shell scripts:

`set -e` and `set -x` only work in shell scripts - they're bash commands. You need to use other approaches in Python and R.

## Automation with Snakemake!

Automation via shell script is wonderful, but there are a few problems here.

First, you have to run the entire workflow each time and it recomputes everything every time. If you're running a workflow that takes 4 days, and you change a command at the end, you'll have to manually go in and just run the stuff that depends on the changed command.

Second, it's very _explicit_ and not very _generalizable_. If you want to run it on a different RNAseq data set, you're going to have to change a lot of commands.

snakemake is one of several workflow systems that help solve these problems. [(You can read the documentation here.)](https://snakemake.readthedocs.io/en/stable/) Let's take a look!

First, let's activate our snakemake environment
```
source deactivate
source activate snake
```

We're going to automate the same script for trimming, but in snakemake. 

```
rule all:
    input:
        "trim/TARA_135_SRF_5-20_rep1_1m_1.qc.fq.gz",
        "trim/TARA_135_SRF_5-20_rep1_1m_2.qc.fq.gz"

rule trim_reads:
    input:
        r1="data/TARA_135_SRF_5-20_rep1_1m_1.fq.gz",
        r2="data/TARA_135_SRF_5-20_rep1_1m_2.fq.gz",
        adapters="trim/combined.fa"
    ouput:
        p1="trim/TARA_135_SRF_5-20_rep1_1m_1.qc.fq.gz",
        p2="trim/TARA_135_SRF_5-20_rep1_1m_2.qc.fq.gz",
        s1="trim/TARA_135_SRF_5-20_rep1_1m_1_s1.qc.fq.gz",
        s2="trim/TARA_135_SRF_5-20_rep1_1m_2_s2.qc.fq.gz"
    shell:'''
    trimmomatic PE {input.r1} \
              {input.r2} \
     {output.p1} {output.s1} \
     {output.p2} {output.s2} \
     ILLUMINACLIP:{input.adapters}:2:40:15 \
     LEADING:2 TRAILING:2 \
     SLIDINGWINDOW:4:2 \
     MINLEN:25
    '''
```

Now we can run this

```
cd $PROJECT
snakemake
```

you should see "Nothing to be done."

That's because the trimmed files already exist!

Let's fix that:

```
rm trim/TARA_135_SRF_5-20_rep1*
```

and now, when you run `snakemake`, you should see the Trimmomatic being run. Yay w00t! Then if you run `snakemake` again, you will see that it doesn't need to do anything - all the files are "up to date".

### Adding and environment

We've been using conda environments throughout our workshop. We showed you have to export your tara environment in the bioconda lesson using
`conda env export -n tara -f $PROJECT/tara_conda_environment.yaml`.
We can use this environment in our snakemake rule as well!


```
rule all:
    input:
        "trim/TARA_135_SRF_5-20_rep1_1m_1.qc.fq.gz",
        "trim/TARA_135_SRF_5-20_rep1_1m_2.qc.fq.gz"

rule trim_reads:
    input:
        r1="data/TARA_135_SRF_5-20_rep1_1m_1.fq.gz",
        r2="data/TARA_135_SRF_5-20_rep1_1m_2.fq.gz",
        adapters="trim/combined.fa"
    ouput:
        p1="trim/TARA_135_SRF_5-20_rep1_1m_1.qc.fq.gz",
        p2="trim/TARA_135_SRF_5-20_rep1_1m_2.qc.fq.gz",
        s1="trim/TARA_135_SRF_5-20_rep1_1m_1_s1.qc.fq.gz",
        s2="trim/TARA_135_SRF_5-20_rep1_1m_2_s2.qc.fq.gz"
    conda: "tara_conda_environment.yaml"
    shell:'''
    trimmomatic PE {input.r1} \
              {input.r2} \
     {output.p1} {output.s1} \
     {output.p2} {output.s2} \
     ILLUMINACLIP:{input.adapters}:2:40:15 \
     LEADING:2 TRAILING:2 \
     SLIDINGWINDOW:4:2 \
     MINLEN:25
     '''
```

We aren't going to run this on the cluster right now, because it requires you to be able to download things and we can't do that on nodo. However, this is the syntax to do this in the future. 

### Other resources

We've covered some basics of snakemake today, but if you want another tutorial, we've add one [here](automation.md). 
