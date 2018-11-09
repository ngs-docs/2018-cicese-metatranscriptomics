# Assembling a Metatranscriptome

Learning objectives:

* What is metatranscriptome assembly?
* How do assemblers work?
* Checking the quality of assembly
* Understanding metatranscriptome assembly


The basic idea with any transcriptome assembly is you feed in your reads and you get out a bunch of *contigs* that represent transcripts, or stretches of RNA present in the reads that don't have any long repeats or much significant polymorphism. You run a transcriptome assembly program using the trimmed reads as input and get out a pile of assembled RNA. 

These contigs will represent transcripts that come from the eukaryotic organisms (poly-A mRNAseq reads) found in each environmental sample.

## Install Megahit

We already installed megahit for you, but [setup](setting-up-tara-environment.md), but here's the installation command for future reference.

```
conda install megahit
```

## Link in the trimmed data

We will be using the same set of TARA oceans mRNAseq reads that we trimmed in the last lesson from [Alberti et al., 2017](https://www.nature.com/articles/sdata201793#t1).

Create a new folder `assembly` to work in 

```
cd $PROJECT
mkdir -p assembly
cd assembly
```

Link the khmer-trimmed data we prepared earlier in the newly created folder:
```
ln -fs ${PROJECT}/khmer_trim/*.khmer.pe.fq.gz .
ls
```
## Run the assembler

Let's run an assembly:

```
time megahit --12 TARA_135_SRF_5-20_rep1_1m.khmer.pe.fq.gz,TARA_135_SRF_5-20_rep2_1m.khmer.pe.fq.gz  --memory 8e9 --num-cpu-threads 2 --out-prefix TARA_135_SRF --out-dir ./TARA_135_SRF_khmer -f
```

This will take about 10 minutes; at the end you should see output like this:

```
--- [STAT] 11733 contigs, total 5202861 bp, min 200 bp, max 4235 bp, avg 443 bp, N50 465 bp
--- [Wed Nov  7 02:13:12 2018] ALL DONE. Time elapsed: 431.097547 seconds ---
```

The output assembly will be `TARA_135_SRF_khmer/TARA_135_SRF.contigs.fa`.


## Looking at the assembly

First, let's copy the assembly into our current directory:

```
cp ./TARA_135_SRF_khmer/*contigs.fa tara135_SRF_megahit.fasta
```

Now, look at the beginning:

```
head tara135_SRF_megahit.fasta 
```

These are the transcripts! Yay!


### What can we do with an assembly?

Why would we do an assembly? What advantages does an assembly have
over the reads? And what can we do with this assembly?

* assembly squashes redundant reads, so the assembly should have
  approximately one sequence per transcript, as opposed to the reads,
  which many have many reads per transcript.
  
* assembled contigs are longer than reads, so it is easier to do gene
  search on them. We'll cover this tomorrow!

* assemblies also have fewer errors than the reads do, so sample
  comparisons and so on may be more accurate. However, assembly also may
  eliminate some of the data if it's really low coverage, and abundance
  information is lost as well.


## Further Reference

There are other asssemblers you can try on your data, some of which are listed
on our [References page](references.md). One that we are excited about is [PLASS](plass-paladin.md), 
which does assembly at the protein level!


