

## Why (or why not) do k-mer trimming?

As we just discussed, raw reads contain errors. To properly compare these samples,
we need to remove these errors. This can also speed up assembly and reduce memory 
requirements (although many assemblers have built in k-mer trimming mechanisms as well).

## Set up workspace and install khmer 

[khmer documentation](http://khmer.readthedocs.io/en/latest)


```
conda install khmer
```

Now, let's create a directory to work in:
```
cd ${PROJECT}
mkdir -p khmer_trim
cd khmer_trim
```

And link in the `qc` trimmed files.
```
ln -s ${PROJECT}/trim/*qc.fq.gz ./
```

To run error trimming, use the khmer script `trim-low-abund.py`:

```
for filename in *_1.qc.fq.gz
do
  #Use the program basename to remove _1.qc.fq.gz to generate the base
  base=$(basename $filename _1.qc.fq.gz)
  echo $base

  #Run khmer trimming
  interleave-reads.py ${base}_1.qc.fq.gz ${base}_2.qc.fq.gz | \
  trim-low-abund.py - -V -Z 10 -C 3 -o - --gzip -M 8e9 | \
  extract-paired-reads.py --gzip -p ${base}.khmer.pe.fq.gz -s ${base}.khmer.se.fq.gz

done
```

## Assess changes in kmer abundance

To see how many k-mers we removed, you can examine the distribution as above,
or use the `unique-kmers.py` script. Let's compare kmers for one sample.

```
unique-kmers.py TARA_135_SRF_5-20_rep1_250k_1.qc.fq.gz TARA_135_SRF_5-20_rep1_250k_2.qc.fq.gz
unique-kmers.py TARA_135_SRF_5-20_rep1_250k.khmer.pe.fq.gz
```  

The first two files (the raw inputs) have a total estimated number of
unique 32-mers of 26760613; the second (trimmed) file has a total
estimated number of unique 32-mers 26659070.  So the trimming removed
approximately 100,000 k-mers.

These numbers are virtually the same BECAUSE WE ARE USING SMALL SUBSET
DATA SETS. For any real data sets, the second number will be MUCH smaller
than the first, indicating that many low-abundance k-mers were removed as
likely errors.

## How does kmer trimming impact sample distances?

We can also test whether we see a difference in our sample distance by running sourmash 
compare again. 

Calculate signatures on your newly trimmed reads:

```
for infile in *pe.fq.gz
do
    sourmash compute -k 31 --scaled 10000 --track-abundance -o ${infile}.sig ${infile}
done
```

Compare these signatures. Let's save the output as a csv to recreate the MDS plot we made before
```
sourmash compare -k 31 --csv tara.comp.csv *sig 
```

And make another MDS plot:

```
wget https://raw.githubusercontent.com/ngs-docs/2018-cicese-metatranscriptomics/master/scripts/mds_plot.R
Rscript mds_plot.R tara.comp.csv tara-comp-mds.pdf 
```

You can look at the output [here](https://github.com/ngs-docs/2018-cicese-metatranscriptomics/blob/master/docs/files/compare-mds-plot.pdf).

Are these samples any different than the untrimmed reads?





For more info and some neat visualization of kmer trimming, go [here](https://2017-cicese-metagenomics.readthedocs.io/en/latest/kmer_trimming.html)
