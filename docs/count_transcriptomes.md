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
