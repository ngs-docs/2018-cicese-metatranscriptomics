# Short Read Quality Control

_lesson changes in progress_
todo:
* change download links
* change names in all downstream steps
* add back in multiqc

You should now be logged into your cluster! If not, go back 
to [Cluster Login](cicese-cluster.md)

## Data source

We will be using RNAseq reads from a small subset of data from the [TARA Oceans Expedition](https://oceans.taraexpeditions.org/en/m/about-tara/les-expeditions/tara-oceans) 
and analyzed as part of [A global ocean atlas of eukaryotic genes](https://www.nature.com/articles/s41467-017-02342-1).

## Set up workspace and download the data 

First, make some directories to work in:
```
cd
mdkir -p tara/data
```

Next, change into the data dir and download the data subsets:
```
cd tara/data
#curl -L https://osf.io/365fg/download -o tara135_5-20_1m.zip
#curl -L https://osf.io/9tf2g/download -o tara136-137_5-20_1m.zip 
#unzip tara135_5-20_1m.zip
#unzip tara136-137_5-20_1m.zip
```

To make our lives easier, define a variable for the location
of this tara working directory:
```
 export PROJECT= ~/tara
```

Check that your data is where it should be

```
ls $PROJECT/data
```

If you see all the files you think you should, good!  Otherwise, debug.

These are FASTQ files -- let's take a look at them:

```
#less 0Hour_ATCACG_L002_R1_001.fastq
```

(use the spacebar to scroll down, and type 'q' to exit 'less')

Question:

* where does the filename come from?
* why are there 1 and 2 in the file names?

Links:

* [FASTQ Format](http://en.wikipedia.org/wiki/FASTQ_format)



## Quality trimming and light quality filtering

Make sure you've got the PROJECT location defined, and your data is there:

```
set -u
printf "\nMy raw data is in $PROJECT/data/, and consists of $(ls -1 ${PROJECT}/data/*.fastq | wc -l) files\n\n"
set +u
```
*Important:* If you get an error above or the count of files is wrong...STOP!! Revisit the installation instructions!

# FROM SEATAC LESSON:
### Link your data into your working directory

Change into your project directory and make a workspace for quality trimming:

```
cd ${PROJECT}
mkdir -p quality
cd quality
```

Now, link the data files into your new workspace

```
ln -s ../data/*.fastq .
```

(Linking with `ln` avoids having to make a copy of the files, which will take up storage space.)

Check to make sure it worked

```
printf "I see $(ls -1 *.fastq | wc -l) files here.\n"
```

You can also do an ``ls`` to list the files.

If you see only one entry, `*.fastq`, then the ln command above didn't work properly.  One possibility is that your files aren't in your data directory; another is that their names don't end with
`.fastq`.

### FastQC


We're going to use
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
summarize the data. We already installed 'fastqc' on our computer for
you.

Now, run FastQC on two files:

```
fastqc *.fastq
```

After this finishes running (has to run on each file so might take a while), type 'ls':

```
ls -d *fastqc.zip*
```

to list the files, and you should see a number of files with the extensions `.fastqc.zip`.

Inside each of the fatqc directories you will find reports from the fastqc. You can download these files using your RStudio Server console, if you like. To install and run an RStudio Server, go [here](https://angus.readthedocs.io/en/2017/visualizing-blast-scores-with-RStudio.html#installing-and-running-rstudio-on-jetstream).

Questions:

* What should you pay attention to in the FastQC report?
* Which is "better", file 1 or file 2? And why?

Links:

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQC tutorial video](http://www.youtube.com/watch?v=bz93ReOv87Y)
* [Examples of fastqc after technical sequencer problems](http://angus.readthedocs.io/en/2015/_static/2015-lecture2-sequencing.pptx.pdf)(starting on slide 40)

There are several caveats about FastQC - the main one is that it only calculates certain statistics (like duplicated sequences) for subsets
of the data (e.g. duplicate sequences are only analyzed for the first 100,000 sequences in each file.

### MultiQC


If you would like to aggregate all of your fastqc reports across many samples, [MultiQC](http://multiqc.info/) will do this into a single report for easy comparison.

Run MultiQC:

```
multiqc .
```

The terminal output should look like this:

```
[INFO   ]         multiqc : This is MultiQC v1.6
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
[INFO   ]          fastqc : Found 20 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```

You can also view this output [here](_static/multiqc_report.html)

#### View your files on your own computer
As an alternative to viewing the files on the Rstudio server, we can secure copy (scp) these files to our own laptops, and view them from there.
```
mkdir ~/Desktop/nema_fastqc  # make a directory for these files
scp username@ip.address:/work/quality/*html ~/Desktop/nema_fastqc
```
where the first argument after `scp` is your login and path for the files we want to copy (from the jetstream instance), and the second argument is the path to place the files on our own computer.

If you are unable to use scp though a terminal output, you can see the fastqc html output [here](_static/6Hour_CGATGT_L002_R2_003.extract_fastqc.html)


### Adapter trim each pair of files

Setup trim directory:
```
cd ..
mkdir trim
cd trim
ln -s ../data/*.fastq .
cat /opt/miniconda3/share/trimmomatic*/adapters/* > combined.fa
```

See excellent paper on trimming from [MacManes 2014](http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full).

Run:

```
for filename in *_R1_*.fastq
do
# first, make the base by removing fastq
  base=$(basename $filename .fastq)
  echo $base

# now, construct the R2 filename by replacing R1 with R2
  baseR2=${base/_R1_/_R2_}
  echo $baseR2

# finally, run Trimmomatic
  trimmomatic PE ${base}.fastq ${baseR2}.fastq \
    ${base}.qc.fq.gz s1_se.gz \
    ${baseR2}.qc.fq.gz s2_se.gz \
    ILLUMINACLIP:combined.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:25

# save the orphans
  zcat s1_se.gz s2_se.gz >> orphans.qc.fq.gz
  rm -f s1_se.gz s2_se.gz
done

```

Now, run fastqc again on trimmed files:
```
fastqc *.qc.fq.gz
multiqc .
```

The paired sequences output by this set of commands will be in the files ending in ``.qc.fq.gz``, with any orphaned sequences all together
in ``orphans.qc.fq.gz``.

Make these trimmed reads read-only and keep them, as we will reuse them later.

```
chmod a-w ${PROJECT}/trim/*.qc.fq.gz
```

Questions:

* How do you figure out what the parameters mean?
* How do you figure out what parameters to use?
* What adapters do you use?
* What version of Trimmomatic are we using here? (And FastQC?)
* Do you think parameters are different for RNAseq and genomic data sets?
* What's with these annoyingly long and complicated filenames?
* why are we running R1 and R2 together?

For a discussion of optimal trimming strategies, see
[MacManes, 2014](http://journal.frontiersin.org/Journal/10.3389/fgene.2014.00013/abstract)
-- it's about RNAseq but similar arguments should apply to metagenome
assembly.

Links:

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

Questions:

* is the quality trimmed data "better" than before?
* Does it matter that you still have adapters!?
