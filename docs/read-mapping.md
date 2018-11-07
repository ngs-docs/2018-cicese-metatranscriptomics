# Mapping reads to a metatranscriptome

We will use [Salmon](http://salmon.readthedocs.org/en/latest/) to
quantify expression. Salmon is a new breed of software for quantifying RNAseq reads that is both really fast and takes
transcript length into consideration ([Patro et al. 2017](https://www.nature.com/articles/nmeth.4197)).

For further reading, see

  * Intro blog post: http://robpatro.com/blog/?p=248
  * A 2016 blog post evaluating and comparing methods [here](https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/)
  * Salmon github repo [here](https://github.com/COMBINE-lab/salmon)
  * 2018 paper: [A direct comparison of genome alignment and transcriptome pseudoalignment](https://www.biorxiv.org/content/early/2018/10/16/444620)


## Quantification with Salmon

Check that salmon is available and see run options:
```
salmon -h
```

Now let's check that we still have the trimmed data we created day 1:

```
set -u
printf "\nMy trimmed data is in $PROJECT/quality/, and consists of $(ls -1 ${PROJECT}/quality/*.qc.fq.gz | wc -l) files\n\n"
set +u

```
where set -u should let you know if you have any unset variables, i.e. if the `$PROJECT` variable is not defined. 

If you see `-bash: PROJECT: unbound variable`, then you need to set the $PROJECT variable.  
```
export PROJECT= ~/work
```
and then re-run the `printf` code block.

NOTE: if you do not have files, please rerun quality trimming steps [here](short-read-quality-control.md)

## Create a directory to work in

```
   cd ${PROJECT}
   mkdir -p quant
   cd quant
```

## Link an assembly

We link a full assembly to use for mapping. This assembly was made with all TARA_135 SRF reads, rather than the subset we used in the [assembly](megahit-assembly.md) tutorial.


```
ln -s /LUSTRE/bioinformatica_data/bioinformatica2018/assembly/tara_f135_full_megahit.fasta ./
```

Note: if you prefer, you can use the annotated assembly we generated with the read subsets instead
```
    #ln -s ${PROJECT}/annotation/tara_f135_SRF.fasta.dammit/tara_f135_SRF.fasta.dammit.fasta ./tara_f135_SRF_annot.fasta 
```

## Run Salmon

Build an index for your new transcriptome:

```
    salmon index --index tara135 --transcripts tara_f135_full_megahit.fasta
```
Next, link in the QC reads (produced in [quality](short-read-quality-trimming.md):

```
   ln -s ../quality/*.qc.fq.gz ./
```


Then, run the salmon command:

```
for sample in *1.qc.fq.gz
do
  base=$(basename $sample _1.fq.gz)
  echo $base
  salmon quant -i tara135 -p 2 -l A -1 ${base}_1.qc.fq.gz -2 ${base}_2.qc.fq.gz -o ${base}_quant
done
```

This will create a bunch of directories named something like
`TARA_135_SRF_5-20_rep1_quant`, containing a bunch of files. Take a
look at what files there are:

You should see:

```
   aux_info
   cmd_info.json
   lib_format_counts.json
   libParams
   logs
   quant.sf
```
The two most interesting files are `quant.sf`, which contains the counts,
and `salmon_quant.log` (in the `logs` directory), which contains a
log from the salmon run.

# Working with the counts

The `quant.sf` files actually contain the relevant information about
expression -- take a look

```
   head -20 TARA_135_SRF_5-20_rep1_1m_quant/quant.sf 
```

You should see output that looks like this:
```
Name	Length	EffectiveLength	TPM	NumReads
k119_5	212	63.757	0.000000	0.000
k119_10	231	75.730	0.000000	0.000
k119_14	261	97.683	0.000000	0.000
k119_16	301	130.690	11.736728	1.000
k119_18	302	131.560	0.000000	0.000
k119_21	203	58.357	0.000000	0.000
k119_22	308	136.790	0.000000	0.000
```



The first column contains the transcript names, and the
fifth column is the "raw counts", which is what many 
differential expression programs need.


# Other useful tutorials and references
https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/salmon.rst
http://angus.readthedocs.io/en/2016/rob_quant/tut.html
https://2016-aug-nonmodel-rnaseq.readthedocs.io/en/latest/quantification.html
