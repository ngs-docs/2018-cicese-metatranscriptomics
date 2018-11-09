
Let's do some basic evaluation of our assembly.


## Set up the directory 

First, make sure you have the `PROJECT` variable defined:

```
echo $PROJECT
```
if you don't see any output, make sure to redefine the $PROJECT variable.

	
Now, let's create a folder `evaluation` to work in:

```
cd $PROJECT
mkdir -p evaluation
cd evaluation
```

Now, let's copy in the assembly we made in the previous lesson:

```
cp $PROJECT/assembly/tara135_SRF_megahit.fasta ./
```

## Generate assembly statistics with Transrate

[Transrate](http://hibberdlab.com/transrate/getting_started.html) is a program that can be used for a couple different types of assembly evaluation. The first and simplest method is to calculate length-based metrics about the assembly, such as the total number of bases, and the N50 of the contigs. Transrate can also be used to compare two assemblies or give you a score which represents proportion of input reads that provide positive support for the assembly. For a further explanation of metrics and how to run the reference-based transrate, see the [documentation](http://hibberdlab.com/transrate/metrics.html) and the paper by [Smith-Unna et al. 2016](http://genome.cshlp.org/content/early/2016/06/01/gr.196469.115). 

We have installed transrate for you. However, if you need to install it in the future, see installation instructions [here](setting-up-tara-environment.md).

See options for running transrate:
```
transrate -h
```

Let's use transrate to calculate some stats on our assembly contigs:

```
transrate --assembly tara135_SRF_megahit.fasta
```

You should see output that looks like this:
```
[ INFO] 2018-11-06 23:50:35 : Loading assembly: /LUSTRE/bioinformatica_data/bioinformatica2018/assembly/tara135_SRF_megahit.fasta
[ INFO] 2018-11-06 23:50:35 : Analysing assembly: /LUSTRE/bioinformatica_data/bioinformatica2018/assembly/tara135_SRF_megahit.fasta
[ INFO] 2018-11-06 23:50:35 : Results will be saved in /LUSTRE/bioinformatica_data/bioinformatica2018/assembly/transrate_results/tara135_SRF_megahit
[ INFO] 2018-11-06 23:50:35 : Calculating contig metrics...
[ INFO] 2018-11-06 23:50:35 : Contig metrics:
[ INFO] 2018-11-06 23:50:35 : -----------------------------------
[ INFO] 2018-11-06 23:50:35 : n seqs                         1502
[ INFO] 2018-11-06 23:50:35 : smallest                        200
[ INFO] 2018-11-06 23:50:35 : largest                        4998
[ INFO] 2018-11-06 23:50:35 : n bases                      638347
[ INFO] 2018-11-06 23:50:35 : mean len                      425.0
[ INFO] 2018-11-06 23:50:35 : n under 200                       0
[ INFO] 2018-11-06 23:50:35 : n over 1k                        40
[ INFO] 2018-11-06 23:50:35 : n over 10k                        0
[ INFO] 2018-11-06 23:50:35 : n with orf                      331
[ INFO] 2018-11-06 23:50:35 : mean orf percent              83.54
[ INFO] 2018-11-06 23:50:35 : n90                             232
[ INFO] 2018-11-06 23:50:35 : n70                             360
[ INFO] 2018-11-06 23:50:35 : n50                             453
[ INFO] 2018-11-06 23:50:35 : n30                             599
[ INFO] 2018-11-06 23:50:35 : n10                             935
[ INFO] 2018-11-06 23:50:35 : gc                             0.51
[ INFO] 2018-11-06 23:50:35 : bases n                           0
[ INFO] 2018-11-06 23:50:35 : proportion n                    0.0
[ INFO] 2018-11-06 23:50:35 : Contig metrics done in 0 seconds
[ INFO] 2018-11-06 23:50:35 : No reads provided, skipping read diagnostics
[ INFO] 2018-11-06 23:50:35 : No reference provided, skipping comparative diagnostics
[ INFO] 2018-11-06 23:50:35 : Writing contig metrics for each contig to /LUSTRE/bioinformatica_data/bioinformatica2018/assembly/transrate_results/tara135_SRF_megahit/contigs.csv
[ INFO] 2018-11-06 23:50:35 : Writing analysis results to assemblies.csv
```

## Comparing Assemblies

We built a metatranscriptome with the full set of TARA_SRF reads. Copy this into your evaluation directory
```
cd ${PROJECT}/evaluation
cp tara125_SRF_full_megahit.fasta ./
```

* How do the two transcriptomes compare with each other?

```
transrate --reference=tara135_SRF_full_megahit.fasta --assembly=tara135_SRF_megahit.fasta --output=full_v_subset
transrate --reference=tara135_SRF_megahit.fasta --assembly=tara135_SRF_full_megahit.fasta --output=subset_v_full
```

## Assess Read Mapping

It's useful to know how well the transcripts represent the sequenced reads. To do this, we'll need to link in the reads we used to generate this assembly: 

```
ln -s ${PROJECT}/trimmed/TARA_135_SRF_5-20_*.qc.fq.gz ./
```

Transrate actually  has a `read assessment` mode that uses `salmon` to "align" reads to the transcriptome and generates some metrics on read mapping. 

We won't run this today, but you could run it via:

```
transrate --assembly=tara135_SRF_megahit.fasta --threads=2 --left=*_1.qc.fq.gz --right *_2.qc.fq.gz --output=${PROJECT}/evaluation/tara135_SRF_transrate
```
