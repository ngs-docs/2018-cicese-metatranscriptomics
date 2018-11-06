# K-mer Error Trimming

(Optional)

khmer documentation: http://khmer.readthedocs.io/en/latest

## Why (or why not) do k-mer trimming?

If you can assemble your data set without k-mer trimming, there's no
reason to do it.  The reason we're error trimming here is to speed up
the assembler (by removing data) and to decrease the memory requirements
of the assembler (by removing a number of k-mers).

## Set up workspace and install khmer 

```
conda install khmer
```

To run error trimming, use the khmer script `trim-low-abund.py`:

```
cd ${PROJECT}/trim

for filename in *_1.qc.fq.gz
do
  #Use the program basename to remove _1.qc.fq.gz to generate the base
  base=$(basename $filename _1.qc.fq.gz)
  echo $base

  #Run khmer trimming
  (interleave-reads.py ${base}_1.qc.fq.gz ${base}_2.qc.fq.gz )| \
  (trim-low-abund.py - -V -Z 10 -C 3 -o - --gzip -M 8e9) | \ 
  (extract-paired-reads.py --gzip -p ${base}.khmer.pe.fq.gz -s ${base}.khmer.se.fq.gz)

done
```

## Assess changes in kmer abundance

To see how many k-mers we removed, you can examine the distribution as above,
or use the `unique-kmers.py` script. Let's compare kmers for one sample.

```
unique-kmers.py TARA_135_SRF_5-20_rep1_1m_1.qc.fq.gz TARA_135_SRF_5-20_rep1_1m_2.qc.fq.gz
unique-kmers.py TARA_135_SRF_5-20_rep1_1m.khmer.pe.fq.gz
```  

