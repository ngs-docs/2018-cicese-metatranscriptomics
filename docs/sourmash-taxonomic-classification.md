# Taxonomic Classification with Sourmash

Knowing what species are in our metatranscriptome allows us to take
advantage of other resources already developed for those species. 
We will use a tool called [sourmash](https://sourmash.readthedocs.io/en/latest/).
Sourmash works by calculating a "signature" for a sequence. For taxonomic 
classification, you can then compare that signature to signatures made
from publicly available data. To make this comparison faster, we put the
publicly available sequences in a database that we can compare against. 

We have made databases from eukaryotic RNA sequences available in GenBank, 
as well as the MMETSP transcriptomes that 
[were recently reassembled](https://figshare.com/articles/Marine_Microbial_Eukaryotic_Transcriptome_Sequencing_Project_re-assemblies/3840153). 
After making signatures for our raw reads, we will search against all of 
our databases at once. 

sourmash cannot find matches across large evolutionary distances.
sourmash seems to work well to search and compare data sets for matches 
at the species and genus level, but does not have much sensitivity beyond 
that. It seems to be particularly good at strain-level analysis. You should 
use protein-based analyses to do searches across larger evolutionary distances
(e.g. with a tool like [KEGG GhostKOALA](https://www.kegg.jp/ghostkoala/)).

We will use the raw reads to do taxonomic classification. Although the errors
present in the raw reads may slow sourmash down a little, there is a chance 
when error trimming transcriptomes that real biological variation is removed. 
Because sourmash is so specific, the errors don't give false positives and so 
keeping all of the possible real variation could improve taxonomic recall in 
some cases.

Let's get started!

First, let's make a directory and link in our raw reads.
```
mkdir -p
ln -s ... .
```

Next, we'll compute a signature of our raw reads. 

sourmash uses k-mers to calculate signatures. A k-mer size of 21 is approximately
specific at the genus level, a 31 is at the species level, and 51 at the strain 
level. We will calculate our signature with all three k-mer sizes so we can 
choose which one we want to use later. 

The `--scaled` parameter in the command tells sourmash to only look at 1/2000th of 
the unique k-mer space. This parameter is agnostic to sequence content, but is 
consistent across signatures. 

Lastly, we use the `--track-abundance` flag. This tells sourmash to keep track of
how often it sees a k-mer. This is meaningful for transcriptomic data because 
expression levels are variable.  

```
sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o ... ...
```

Now let's download some databases that contain signatures from transcriptomes of 
publicly-available eukaryotic sequences. We will download databases for plants, 
fungi, protozoa, invertebrates, and vertebrates, as well a database that contains 
sequences from the MMETSP re-assembly project (>700 marine eukaryotic transcriptomes). 
Even if you don't think a certain class of eukaryote is likely to be in your sample
(for instance, vertebrate_mammalian), it's a good idea to include all of the databases
anyway. It can help you quickly locate potential contamination in your sequencing (we
wouldn't expect to find mouse RNA in the ocean, but it is a commonly sequenced organism
and thus if you find it, if may be a contaminant). 

```
mkdir -p sourmash_databases
cd sourmash_databases
wget -O genbank-rna-vertebrate_other-k31.tar.gz https://osf.io/qgyax/download
wget -O genbank-rna-vertebrate_mammalian-k31.tar.gz https://osf.io/6c9uy/download
wget -O genbank-rna-invertebrate-k31.tar.gz https://osf.io/7v8ck/download
wget -O ...
wget -O ...
wget -O ...
```


And now we can perform taxonomic classification!

```
sourmash gather -k 31 --scaled 2000 sig db(s)
```


