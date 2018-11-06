# Short Read Quality Control

still to do:  
  * viewing fastqc files (use Rstudio server to view files still? OR scp?)

You should now be logged into your cluster! If not, go back 
to [Cluster Login](cicese-cluster.md)

## Data source

We will be using RNAseq reads from a small subset of data from the [TARA Oceans Expedition](https://oceans.taraexpeditions.org/en/m/about-tara/les-expeditions/tara-oceans), 
from [Alberti et al., 2017](https://www.nature.com/articles/sdata201793#t1) and analyzed as part of [A global ocean atlas of eukaryotic genes](https://www.nature.com/articles/s41467-017-02342-1).

## Set up workspace and download the data 

First, make some directories to work in:
```
cd
mdkir -p work/data
```


Next, change into the data dir and download the data subsets:
```
cd work/data
curl -L https://osf.io/76s3r/download -o tara135_5-20_1m.zip
curl -L https://osf.io/smcyf/download -o tara136-137_5-20_1m.zip 
```

Now, let's unzip and make the files difficult to delete
```
unzip tara135_5-20_1m.zip
unzip tara136-137_5-20_1m.zip

chmod u-w */*fq.gz
```

To make life easier, let's define a variable for the location
of this tara working directory:
```
 export PROJECT=~/work
```

Check that your data is where it should be

```
ls $PROJECT/data/*/
```

If you see all the files you think you should, good!  Otherwise, debug.

These are FASTQ files -- let's take a look at them:

```
zless $PROJECT/data/tara135_5-20_1m/TARA_135_DCM_5-20_rep1_1m_1.fq.gz
```
(use the spacebar to scroll down, and type 'q' to exit 'zless')

Question:

* where does the filename come from?
* why are there 1 and 2 in the file names?

Links:

* [FASTQ Format](http://en.wikipedia.org/wiki/FASTQ_format)



## Quality trimming and light quality filtering

Make sure you've got the PROJECT location defined, and your data is there:

```
set -u
printf "\nMy raw data is in $PROJECT/data/, and consists of $(ls -1 ${PROJECT}/data/*/*.fq.gz | wc -l) files\n\n"
set +u
```
*Important:* If you get an error above or the count of files is wrong...STOP!! Revisit the download & unzip instructions!

### Link your data into your working directory

Change into your project directory and make a workspace for quality trimming:

```
cd ${PROJECT}
mkdir -p quality
cd quality
```

Now, link the data files into your new workspace

```
ln -s ../data/*/*.fq.gz ./
```

(Linking with `ln` avoids having to make a copy of the files, which will take up storage space.)

Check to make sure it worked

```
printf "I see $(ls -1 *.fq.gz | wc -l) files here.\n"
```

You can also do an ``ls`` to list the files.

If you see only one entry, `*.fq.gz`, then the ln command above didn't work properly.  One possibility is that your files aren't in your data directory; another is that their names don't end with
`.fq.gz`.

### FastQC


We're going to use
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to summarize the data. 

To install fastqc via conda:
Note: make sure you've followed the [conda setup instructions](working-with-bioconda.md).
```
conda install fastqc
```

Now, run FastQC on the files:

```
fastqc *.fq.gz
```

After this finishes running (has to run on each file so might take a while), type 'ls':

```
ls -d *fastqc.zip*
```

to list the files, and you should see a number of files with the extensions `.fastqc.zip`.

Inside each of the fastqc directories you will find reports from the fastqc. You can download these files using your RStudio Server console, if you like. To install and run an RStudio Server, go [here](https://angus.readthedocs.io/en/2017/visualizing-blast-scores-with-RStudio.html#installing-and-running-rstudio-on-jetstream).

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

Install MultiQC with conda:

```
conda install multiqc
```

Run MultiQC:

```
multiqc .
```

The terminal output should look like this:

```
[INFO   ]         multiqc : This is MultiQC v1.6
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
[INFO   ]          fastqc : Found 16 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```

You can also view this output [here](files/multiqc_report.html)

#### View your files on your own computer
As an alternative to viewing the files on the Rstudio server, we can secure copy (scp) these files to our own laptops, and view them from there.
```
mkdir ~/Desktop/tara_fastqc  # make a directory for these files
scp username@ip.address:~/work/quality/*html ~/Desktop/tara_fastqc
```
where the first argument after `scp` is your login and path for the files we want to copy (from the jetstream instance), and the second argument is the path to place the files on our own computer.

If you are unable to use scp though a terminal output, you can see the fastqc html output [here](files/TARA_135_DCM_5-20_rep1_1m_1_fastqc.html)


### Adapter trim each pair of files

Install Trimmomatic:
```
conda install trimmomatic
```

Setup trim directory:
```
cd $PROJECT
mkdir -p trim
cd trim
ln -s ../data/*/*.fq.gz .
cat ~/miniconda3/envs/tara/share/trimmomatic*/adapters/* > combined.fa
```

See excellent paper on trimming from [MacManes 2014](http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full).

Run:

```
for filename in *1.fq.gz
do
#Use the program basename to remove _1.fq.gz to generate the base
base=$(basename $filename _1.fq.gz)
echo $base

# run Trimmomatic
trimmomatic PE ${base}_1.fq.gz \
              ${base}_2.fq.gz \
     ${base}_1.qc.fq.gz ${base}_s1_se \
     ${base}_2.qc.fq.gz ${base}_s2_se \
     ILLUMINACLIP:combined.fa:2:40:15 \
     LEADING:2 TRAILING:2 \
     SLIDINGWINDOW:4:2 \
     MINLEN:25
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
