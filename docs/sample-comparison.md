# Comparing Metatranscriptome samples

Often times, the goal of a (meta)transcriptome sequencing project is
to compare the differences in functional profiles between samples. One way
to do this is with differential expression. To do differential expression,
we de novo assemble a metatranscriptome, quantify the number of reads that "map"
back to the metatranscriptome, and use those counts to find differences. However,
assembly can be a long process, and often we want to get to know our data a bit
before we launch into that process. 

We can use k-mer profiles of the reads to compare samples using `sourmash compare`.

First, let's make a directory that we will be working in:
```
mkdir -p ${PROJECT}/sourmash-compare
cd ${PROJECT}/sourmash-compare
```

Let's link in our trimmed reads. 

```
ln -s $PROJECT/trim/*1.qc.fq.gz .
ln -s $PROJECT/trim/*2.qc.fq.gz . 
```

Now we can calculate signatures for each of the files. This will take 5 or 10 minutes to run

```
for infile in *1.qc.fq.gz
do
    j=$(basename ${infile} 1.qc.fq.gz)
    sourmash compute -k 31 --scaled 10000 --track-abundance --merge ${j} -o ${j}.sig ${j}2.qc.fq.gz ${infile}
done
```

Using these signatures, we can compare our samples. 

```
sourmash compare -k 31 -o tara-trimmed.comp *sig
```

Now let's plot! Sourmash has a built in plot utility that we can take advantage of.
The output is a heatmap. 

```
sourmash plot --labels tara-trimmed.comp
```

This produces a file, `tara-trimmed.comp.matrix.png`, that contains a
similarity matrix.
You can see the heatmap [here](https://github.com/ngs-docs/2018-cicese-metatranscriptomics/blob/master/docs/files/tara-trimmed.comp.matrix.png)

We can also use the output of sourmash compare to calculate an MDS plot. Let's 
rerun `sourmash compare`, this time saving the output in csv format.
```
sourmash compare -k 31 --csv tara-trimmed.comp.csv *sig 
```

We can use this output to make a Multidimensional Scaling plot. MDS plots are 
commonly used in transcriptome workflows to visualize distance between samples. 
Here the strength is we used all of our reads to calculate these distances. 

To make an MDS plot, run:
```
ln -s /LUSTRE/bioinformatica_data/bioinformatica2018/scripts/mds_plot.R . 
Rscript mds_plot.R tara-trimmed.comp.csv tara-trimmed-comp-mds.pdf 
```

The script source is [here](https://raw.githubusercontent.com/ngs-docs/2018-cicese-metatranscriptomics/master/scripts/mds_plot.R) if you are interested!

This outputs a file `tara-trimmed-comp-mds.pdf`.  You can see what
that visualization looks like
[here](https://github.com/ngs-docs/2018-cicese-metatranscriptomics/blob/master/docs/files/tara-trimmed-comp-mds.pdf).

We see that our samples cluster by site (`TARA_135` vs `TARA_136`) and
then by depth (SRF for surface vs DCM for deep cholorophyl maximum).

Throughout this lesson we have been working with raw reads. 
Raw reads contain a lot of errors, and these errors are included in the
signatures. Next we will learn to k-mer trim our reads, and then re-run 
`compare` to see if makes a difference!
