# Assembling with PLASS


The basic idea with any transcriptome assembly is you feed in your reads and you get out a bunch of *contigs* that represent transcripts, or stretches of RNA present in the reads that don't have any long repeats or much significant polymorphism. You run a transcriptome assembly program using the trimmed reads as input and get out a pile of assembled RNA. We used Megahit earlier to build this sort of assembly.

[PLASS](https://plass.mmseqs.org) is a new type of assembler tha works in protein space. First, 

These contigs will represent protein sequences that come from the eukaryotic organisms found in each environmental sample.

## Install PLASS

We already installed plass for you, but [setup](setting-up-tara-environment.md), but here's the installation command for future reference.

```
conda install plass
```

## Link in the trimmed data

We will be using the same set of TARA oceans mRNAseq reads that we trimmed in the last lesson from [Alberti et al., 2017](https://www.nature.com/articles/sdata201793#t1).

Create a new folder `assembly` to work in 

```
cd $PROJECT
mkdir -p assembly_plass
cd assembly_plass
```

Link the khmer-trimmed data we prepared earlier in the newly created folder:
```
ln -fs ${PROJECT}/trim/*.khmer.fq.gz .
ls
```

Plass needs separate files for `_1` and `_2`, so let's split the files

```
for filename in *.khmer.fq.gz
do
  #Use the program basename to remove _1.qc.fq.gz to generate the base
  base=$(basename $filename .khmer.fq.gz)
  echo $base

  #Run khmer trimming
  split-paired-reads.py --gzip ${base}.khmer.fq.gz -1 ${base}_1.khmer.fq.gz -2 ${base}_2.khmer.fq.gz

done
```

## Run the assembler

Let's run an assembly:

```
plass assemble TARA_135_SRF_5-20_rep1_1m_1.khmer.fq.gz TARA_135_SRF_5-20_rep1_1m_2.khmer.fq.gz tara135_srf_plass.fasta --threads 2
```

This will take about XX minutes; at the end you should see output like this:

```
<put some output here>
```

The output assembly will be `tara135_srf_plass.fasta`.

## Looking at the assembly

Let's look at the beginning

```
head tara135_srf_plass.fasta 
```

How is this assembly different from the `megahit` transcripts?


## Mapping to the PLASS Assembly

### Set up the workspace

We'll be using `Paladin` to map back to the plass assembly

Let's make a directory to work in
```
cd $PROJECT
mkdir -p paladin_mapping
cd paladin_mapping
```

Link in the `qc` trimmed reads

```
ln -s ${PROJECT}/trim/*qc.fq.gz ./
```

### Index the Assembly:

```
paladin index -r3  tara135_srf_plass.fasta
```

### Run the mapping:

```
 paladin align -f 125 -t 2 tara135_srf_plass.fasta TARA_135_SRF_5-20_rep1_1m_1.khmer.fq.gz 
```
