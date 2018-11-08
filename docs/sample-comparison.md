# Comparing Metatranscriptome samples

Often times, the goal of a (meta)transcriptome sequencing project is
to compare the differences in functional profiles between samples. One way
to do this is with differential expression. To do differential expression,
we de novo assemble a transcriptome, quantify the number of reads that "map"
back to the transcriptome, and use those counts to find differences. However,
as we have seen throughout this workshop, not even the majority of reads in 
our samples assemble. How can we compare expression profiles when we are only
looking at ~30% of the reads?

We can use k-mer profiles of the reads to compare samples using `sourmash compare`.
During the `gather` lesson, we showed you how to calculate signatures for an assembly.
The process is the same for read files. However, unlike when we ran `gather` with raw
reads, we want to remove k-mers from our samples that are most likely errors before
we compare our samples to each other (with gather, if a k-mer is an error, it is
marked as uncharacterized; with compare, errors could lead samples to being less
similar than they really are). 

Because we already showed you how to error trim reads during the quality control 
section, we have made pre-trimmed and subsampled files available for you to download.
We will use these to calculate signatures with. 

First, let's make a directory and download the k-mer trimmed reads for each of our samples:
```
mkdir -p ${PROJECT}/sourmash-compare
cd ${PROJECT}/sourmash-compare
wget -O tara135.zip https://osf.io/zvnug/download
wget -O tara136-137.zip https://osf.io/7eydv/download
```

Our read files are instead of the zipped files, so we need to uncompress them
```
for infile in *zip
do
   unzip $infile
done
```

Now we can calculate signatures for each of the files. This will take 5 or 10 minutes to run

```
for infile in *.fq.gz
do
    sourmash compute -k 31 --scaled 10000 --track-abundance -o ${infile}.sig ${infile}
done
```

We can now compare our signatures

```
sourmash compare -k 31 -o tara.comp *sig
```

Now let's plot! Sourmash has a built in plot utility that we can take advantage of.
The output is a heatmap. 

```
sourmash plot --labels tara.comp
```

We can also use the output of sourmash compare to calculate an MDS plot. Let's 
rerun `sourmash compare`, this time saving the output in csv format.
```
sourmash compare -k 31 --csv tara.comp.csv *sig 
```

We can use this output to make a Multidimensional Scaling plot. MDS plots are 
commonly used in transcriptome workflows to visualize distance between samples. 
Here the strength is we used all of our reads to calculate these distances. 

To make an MDS plot, run:
```
wget https://raw.githubusercontent.com/ngs-docs/2018-cicese-metatranscriptomics/master/scripts/mds_plot.R
Rscript mds_plot.R tara.comp.csv tara-comp-mds.pdf 
```

You can see what that visualization looks like [here](https://github.com/ngs-docs/2018-cicese-metatranscriptomics/blob/master/docs/files/compare-mds-plot.pdf).

We see that our samples cluster by site and then by depth. We will see if differential
expression recapitulates this pattern!
