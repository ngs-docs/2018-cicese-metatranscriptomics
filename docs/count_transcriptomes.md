In metagenomics, we use binning to approximate the number of genomes (or species) in a sample.
Although we can apply taxonomic classification methods like sourmash to metatranscriptomic 
samples, these methods only work reliably when sequences already exist for the organism of 
interest. This is still a relatively rare occurrence in environmental metatranscriptomics. 

We will make use of other data in our samples to approximate the number of species present.
Although our samples were poly-A selected, inevitably a some ribosomal RNA is always sequenced
given its abundance in samples. We can estimate the number of species in our metatranscriptome
by counting the unique number of ribosomal protein sequences we see. There are many ribosomal
proteins one could use for this analysis, and we provide a table of these at the bottom of 
this lesson. We will use three sequences here, and average the number of unique sequences we 
detect between the three sequences to estimate the number of transcriptomes in our sample.


## Other ribosomal proteins

| name | COG     | PFAM    |
|------|---------|---------|
| rpsG | COG0049 | PF00177 |
| RpsB | COG0052 | PF00318 |
|      | COG0080 |         |
|      | COG0081 |         |
|      | COG0087 |         |
|      | COG0088 |         |
|      | COG0091 |         |
|      | COG0092 |         |
|      | COG0093 |         |
|      | COG0094 |         |
|      | COG0096 |         |
|      | COG0097 |         |
|      | COG0098 |         |
|      | COG0099 |         |
|      | COG0100 |         |
|      | COG0102 |         |
|      | COG0103 |         |
|      | COG0184 |         |
|      | COG0185 |         |
|      | COG0186 |         |
|      | COG0187 |         |
|      | COG0188 |         |
|      | COG0189 |         |
|      | COG0190 |         |
|      | COG0191 |         |
|      | COG0192 |         |
|      | COG0193 |         |
|      | COG0194 |         |
|      | COG0200 |         |
|      | COG0256 |         |
|      | COG0522 |         | 
