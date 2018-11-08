# Comparing Metatranscriptome samples

Often times, the goal of a (meta)transcriptome sequencing project is
to compare the differences in functional profiles between samples. One way
to do this is with differential expression. To do differential expression,
we de novo assemble a transcriptome, quantify the number of reads that "map"
back to the transcriptome, and use those counts to find differences. However,
assembly can be a long process, and often we want to get to know our data a bit
before we launch into that process. 

We can use k-mer profiles of the reads to compare samples using `sourmash compare`.
During the `gather` lesson, we showed you how to calculate signatures for an assembly.
The process is the same for read files. 

First, let's make a directory and download the reads for our samples:
```
mkdir -p ${PROJECT}/sourmash-compare
cd ${PROJECT}/sourmash-compare

# link interleaved reads 
```

Our read files are inside of the zipped files, so we need to uncompress them
```
unzip # add file name
```

Now we can calculate signatures for each of the files. This will take 5 or 10 minutes to run

```
for infile in *.fq.gz
do
    sourmash compute -k 31 --scaled 10000 --track-abundance -o ${infile}.sig ${infile}
done
```

Using these signatures, we can compare our samples. 

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

You can see what that visualization looks like [here]

We see that our samples cluster by site and then by depth. 

However, throughout this lesson we have been working with raw reads. 
Raw reads contain a lot of errors, and these errors are included in the
signatures. Next we will learn to k-mer trim our reads, and then re-run 
`compare` to see if makes a difference!
