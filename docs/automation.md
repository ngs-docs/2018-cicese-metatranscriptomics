# Automation via shell scripts and snakemake

Note: This lesson was written as part of [DIBSI 2018](http://ivory.idyll.org/dibsi/) 
by C. Titus Brown. Links to the data are all provided, so you should be able to follow along. 

 
Learning objectives:

* understand how to write shell scripts
* understand the differences between shell scripts and workflow engines
* develop a first snakemake workflow.

We're going to reprise the RNAseq differential expression, but we're going to automate the heck out of it.

### Clean off old stuff

```
cd ~/
rm -fr data rnaseq
```

### Install the necessary software and data

You'll need to have conda and RStudio installed (see [installation instructions](https://angus.readthedocs.io/en/2018/jetstream-bioconda-config.html)).

Then install salmon and edgeR:

```
cd ~/

conda install -y salmon

curl -L -O https://raw.githubusercontent.com/ngs-docs/angus/2018/scripts/install-edgeR.R
sudo Rscript --no-save install-edgeR.R
```

You'll also need the yeast data from [trimming](https://dibsi-rnaseq.readthedocs.io/en/latest/quality-trimming.html) again:

```
cd ~/
mkdir -p data
cd data

curl -L https://osf.io/5daup/download -o ERR458493.fastq.gz
curl -L https://osf.io/8rvh5/download -o ERR458494.fastq.gz
curl -L https://osf.io/2wvn3/download -o ERR458495.fastq.gz
curl -L https://osf.io/xju4a/download -o ERR458500.fastq.gz
curl -L https://osf.io/nmqe6/download -o ERR458501.fastq.gz
curl -L https://osf.io/qfsze/download -o ERR458502.fastq.gz
```

And, finally, download the necessary analysis code into your home directory

```
cd ~/
curl -L -O https://raw.githubusercontent.com/ngs-docs/angus/2018/scripts/gather-counts.py

curl -L -O https://raw.githubusercontent.com/ngs-docs/angus/2018/scripts/yeast.salmon.R
```

## Write a shell script to run RNAseq analysis

A **shell script** is a text file full of shell commands, that run just as if you're running them interactively at the command line.

Here, we'll create one that takes the commands from the RNAseq differential expression tutorial and runs 'em all at once.

### Put all of the RNAseq commands into a single text file

We'll script everything but the install and data download part of [the RNAseq differential expression tutorial](https://angus.readthedocs.io/en/2018/rna-seq.html).

Use RStudio to create a new file named `run-rnaseq.sh` in your home directory on the Jetstream computer, and write the following commands into it:

```
# make an rnaseq directory and change into it
mkdir rnaseq
cd rnaseq

# link the data in from the other directory
ln -fs ~/data/*.fastq.gz .

# download the reference transcriptome
curl -O https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz

# index the reference transcriptome
salmon index --index yeast_orfs --type quasi --transcripts orf_coding.fasta.gz

# run salmon quant on each of the input read data sets
for i in *.fastq.gz
do
   salmon quant -i yeast_orfs --libType U -r $i -o $i.quant --seqBias --gcBias
done

# collect the counts
python2 ~/gather-counts.py

# run the edgeR script!
Rscript --no-save ~/yeast.salmon.R
```

This is now a shell script that you can use to execute all
of those commands in *one* go -- try it out! Run:

```
cd ~/
bash run-rnaseq.sh
```

You should now have a directory `rnaseq/` containing a bunch of files, including the PDF outputs of the RNAseq pipeline, `yeast-edgeR-MA-plot.pdf`, `yeast-edgeR-MDS.pdf`, and `yeast-edgeR.csv`.

Pretty cool, eh?

### Re-running the shell script

Suppose you wanted to re-run the script. How would you do that?

Well, note that the `rnaseq` directory is created at the top of the script, and everything is executed in that directory. So if you remove the rnaseq directory like so,

```
rm -fr rnaseq
```

you can then do

```
bash run-rnaseq.sh
```

to just ...rerun everything.

### A trick: make it executable

You can get rid of the `bash` part of the command above with
some magic:

Put 
```
#! /bin/bash
```
at the top of the file, and then run

```
chmod +x ~/run-rnaseq.sh
```

at the command line.

You can now run

```
./run-rnaseq.sh
```
instead of `bash run-rnaseq.sh`.

You might be thinking, ok, why is this important? Well, you can do the same with R scripts and Python scripts (but put `/usr/bin/env Rscript` or `/usr/bin/env python` at the top, instead of `/bin/bash`). This basically annotates the script with the language it's written in, so you don't have to know or remember yourself.

So: it's not necessary but it's a nice trick.

You can also always *force* a script to be run in a particular language by specifying `bash <scriptname>` or `Rscript <Scriptname>`, too.

### Another good addition: make it fail nicely

One sort of weird aspect of shell scripts is that (by default) they keep running even if there's an error.  This is bad behavior and we should turn it off.

You can observe this behavior by rerunning the script above without deleting the directory `rnaseq/` - the `mkdir` command will print an error because the directory still exists, but 

A good addition to every shell script is to make it fail on the first error. Do this by putting
```
set -e
``` 
at the top - this tells bash to exit at the first error, rather than continuing bravely onwards.

### A final good addition: make shell scripts print out the commands they're running!

You might notice that the shell script gives you the output of the commands it's running, but doesn't tell you *which* commands it's running.

If you add
```
set -x
```
at the top of the shell script and then re-run it,

```
cd ~/
rm -fr rnaseq
./run-rnaseq.sh
```

then you'll see the full set of commands being run!

### Tracking output of shell scripts

you can do:

```
exec > output.log
exec 2> error.log
```
to save all output from that point on to those files.

Or, you can do:

```
nohup ./run-rnaseq.sh > output.log 2> error.log &
```
and that will do the same thing, while putting it in the background.

You can also run 'screen' or 'tmux' (google it).


Torsten recommends:
```
script session.log
```

While you're at it, you should look up nohup.

### A final note on shell scripts:

`set -e` and `set -x` only work in shell scripts - they're bash commands. You need to use other approaches in Python and R.

## snakemake for automation

Automation via shell script is wonderful, but there are a few problems here.

First, you have to run the entire workflow each time and it recomputes everything every time. If you're running a workflow that takes 4 days, and you change a command at the end, you'll have to manually go in and just run the stuff that depends on the changed command.

Second, it's very _explicit_ and not very _generalizable_. If you want to run it on a different RNAseq data set, you're going to have to change a lot of commands.

snakemake is one of several workflow systems that help solve these problems. [(You can read the documentation here.)](https://snakemake.readthedocs.io/en/stable/) Let's take a look!

### First, install snakemake.

Run:
```
conda install -y snakemake
```

## Writing your first snakemake rule

We're going to automate the R script at the end of the shell
script in snakemake.

So, run the rnaseq pipeline:

```
cd ~/
rm -fr rnaseq
./run-rnaseq.sh
```

Now, using RStudio, create a new text file, and put the following text in it:

```
rule make_plots:
  input:
  output:
    "yeast-edgeR-MA-plot.pdf",
    "yeast-edgeR-MDS.pdf",
    "yeast-edgeR.csv"
  shell:
    "Rscript --no-save ~/yeast.salmon.R"
```
This rule says:
* to make the list of outputs (the PDFs and CSV file),
* you need the list of inputs (currently lbank)
* and then you run the shell command.

Save this with the name `Snakefile` in the `rnaseq/` directory, and run `snakemake make_plots`:

```
cd ~/rnaseq
snakemake make_plots
```

you should see "Nothing to be done."

That's because the PDFs and CSV file already exist!

Let's fix that:

```
rm yeast-*.pdf yeast-*.csv
```

and now, when you run `snakemake make_plots`, you should see the R script being run. Yay w00t! Then if you run `snakemake make_plots` again, you will see that it doesn't need to do anything - all the files are "up to date".



#### Challenge: add the explicit list of inputs.

What does the `make_plots` rule take as input? List them!

Hint: what does the `yeast.salmon.R` script take in?

----

Once you're done, you can also update the timestamp on the counts files and it will regenerate the PDFs --

```
touch *.counts
snakemake make_plots
```
because it knows that the files the rule depends on have changed!

#### Challenge: write a rule that generates the .counts files.

Looking at the shell script, write a new Snakemake rule that generates the counts files as outputs.

You can start by adding a new blank rule:

```
rule make_counts:
  input:
  output:
  shell:
  	""
```

and filling in from there. Note that rule order does not matter in the Snakefile!

Hints:

* you should list each .counts file separately.

To test the rule, do

```
rm *.counts yeast-*.pdf
```

and it will regenerate everything needed to build the PDF!

### Writing rules with wildcards

Let's suppose we want to write a rule that takes a FASTQ file (say, `ERR458501.fastq.gz`) and runs salmon on it to generate the directory `ERR458501.fastq.gz.quant`. How would we do that?

We could write the rule as follows:

```
rule generate_quant_dir_ERR458501:
	input:
		"ERR458501.fastq.gz"
	output:
		directory("ERR458501.fastq.gz.quant")
	shell:
		"salmon quant -i yeast_orfs --libType U -r ERR458501.fastq.gz -o ERR458501.fastq.gz.quant --seqBias --gcBias"
```

and that will work just fine -- you can run,

```
rm -fr ERR458501.fastq.gz.quant
snakemake ERR458501.fastq.gz.quant
```
and it will happily run salmon for you.

But there's a lot of repeated stuff in there. Can we do cool computer things to fix that? Yes! We can use wildcards!

Try replacing `ERR458501` with `{accession}` - this makes it into a _wildcard_ that snakemake can automatically fill in:

```
rule generate_quant_dir:
  input:
    "{accession}.fastq.gz"
  output:
    directory("{accession}.fastq.gz.quant")
  shell:
    "salmon quant -i yeast_orfs --libType U -r {wildcards.accession}.fastq.gz -o {wildcards.accession}.fastq.gz.quant --seqBias --gcBias"
```

Now, you can run

```
snakemake ERR458501.fastq.gz.quant
```
and it will work - but you can _also_ run

```
snakemake ERR458493.fastq.gz.quant
```
and it will do the right thing there, as well.

### Using 'expand' to work with lists of files

The challenge to write a rule that generates the .counts files should have left you with a list of six input files (`.quant` directories), and `make_plots` has a very similar of six input files (`.quant.counts`).

Can we make this more succinct? Yep. You can have a list that you **expand** in various ways.  First, define the list of files at the top of the Snakefile:

```
input_files=["ERR458493.fastq.gz",
  "ERR458494.fastq.gz",
  "ERR458495.fastq.gz",
  "ERR458500.fastq.gz",
  "ERR458501.fastq.gz",
  "ERR458502.fastq.gz"]
```

now, modify the `make_plots` rule to have the input be:

```
rule make_plots:
   input:
      expand("{input_files}.quant.counts", input_files=input_files)
```

#### Challenge:

Update the other rules that have lists of six files using `expand`.

### Visualizing workflow diagrams

Try running:

```
snakemake --dag yeast-edgeR.csv |  dot -Tpng > dag.png
```

### Challenge: Extend the snakemake rules to download and prepare the reference

'nuff said'

### Other things snakemake can do

You can specify what programs your pipeline depends in on your snakemake setup! snakemake interfaces well with conda (in fact, the authors of snakemake are also leads on the bioconda project) so you can pretty much just list your packages.

## Reproducibility and provenance

If you write all your workflows like this, all you need to do is give readers your data files, your Snakefile, and your list of software... and they can reproduce everything you did!

### Answer to first challenge:

```
rule make_plots:
  input:
    "ERR458493.fastq.gz.quant.counts",
    "ERR458500.fastq.gz.quant.counts",
    "ERR458494.fastq.gz.quant.counts", 
    "ERR458501.fastq.gz.quant.counts",
    "ERR458495.fastq.gz.quant.counts",
    "ERR458502.fastq.gz.quant.counts"
  output:
    "yeast-edgeR-MA-plot.pdf",
    "yeast-edgeR-MDS.pdf",
    "yeast-edgeR.csv"
  shell:
    "Rscript --no-save ~/yeast.salmon.R"
```

### Answer to second challenge

```
# make the .counts files
rule make_counts:
  input:
    "ERR458493.fastq.gz.quant",
    "ERR458500.fastq.gz.quant",
    "ERR458494.fastq.gz.quant", 
    "ERR458501.fastq.gz.quant",
    "ERR458495.fastq.gz.quant",
    "ERR458502.fastq.gz.quant"
  output:
    "ERR458493.fastq.gz.quant.counts",
    "ERR458500.fastq.gz.quant.counts",
    "ERR458494.fastq.gz.quant.counts", 
    "ERR458501.fastq.gz.quant.counts",
    "ERR458495.fastq.gz.quant.counts",
    "ERR458502.fastq.gz.quant.counts"
  shell:
    "python2 ~/gather-counts.py"
```

### Final snakemake file

```
input_files=["ERR458493.fastq.gz",
  "ERR458494.fastq.gz",
  "ERR458495.fastq.gz",
  "ERR458500.fastq.gz",
  "ERR458501.fastq.gz",
  "ERR458502.fastq.gz"]
  
rule make_plots:
  input:
    expand("{input_files}.quant.counts", input_files=input_files)
  output:
    "yeast-edgeR-MA-plot.pdf",
    "yeast-edgeR-MDS.pdf",
    "yeast-edgeR.csv"
  shell:
    "Rscript --no-save ~/yeast.salmon.R"
    
    
# make the .counts files
rule make_counts:
  input:
    expand("{input_files}.quant", input_files=input_files)
  output:
    expand("{input_files}.quant.counts", input_files=input_files)
  shell:
    "python2 ~/gather-counts.py"
    
rule generate_quant_dir:
  input:
    "{accession}.fastq.gz"
  output:
    "{accession}.fastq.gz.quant"
  shell:
    "salmon quant -i yeast_orfs --libType U -r {wildcards.accession}.fastq.gz -o {wildcards.accession}.fastq.gz.quant --seqBias --gcBias"
```
