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
After making signatures for our assembly, we will search against all of 
our databases at once. 

sourmash cannot find matches across large evolutionary distances.
sourmash seems to work well to search and compare data sets for matches 
at the species and genus level, but does not have much sensitivity beyond 
that. It seems to be particularly good at strain-level analysis. You should 
use protein-based analyses to do searches across larger evolutionary distances
(e.g. with a tool like [KEGG GhostKOALA](https://www.kegg.jp/ghostkoala/)).


Let's get started!

First, let's make a directory and link in our raw reads.
```
mkdir -p ~/sourmash-gather
cd ~/sourmash-gather
ln -s ... .
```

Next, let's make a signature of our assembly.

sourmash uses k-mers to calculate signatures. A k-mer size of 21 is approximately
specific at the genus level, a 31 is at the species level, and 51 at the strain 
level. We will calculate our signature with all three k-mer sizes so we can 
choose which one we want to use later. 

The `--scaled` parameter in the command tells sourmash to only look at 1/10,000th of 
the unique k-mer space. This parameter is agnostic to sequence content, but is 
consistent across signatures. 

Lastly, we use the `--track-abundance` flag. This tells sourmash to keep track of
how often it sees a k-mer. This is meaningful for transcriptomic data because 
expression levels are variable.  

```
sourmash compute -k 21,31,51 --scaled 10000 --track-abundance -o tara135_SRF_megahit.sig tara135_SRF_megahit.fasta
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
wget -O genbank-rna-vertebrate_other-k31.tar.gz https://osf.io/qgyax/download
wget -O genbank-rna-vertebrate_mammalian-k31.tar.gz https://osf.io/6c9uy/download
wget -O genbank-rna-invertebrate-k31.tar.gz https://osf.io/7v8ck/download
wget -O genbank-rna-fungi-k31.tar.gz https://osf.io/g6mcr/download
wget -O genbank-rna-plant-k31.tar.gz https://osf.io/kctus/download
wget -O genbank-rna-protozoa-k31.tar.gz https://osf.io/fnu2q/download
wget -O mmetsp-k31-named.tar.gz https://osf.io/cdvqn/download
```

Next, we need to uncompress the databases:

```
for infile in *tar.gz
do
    tar xf ${infile}
done
```

And now we can perform taxonomic classification!

```
sourmash gather -k 31 --scaled 10000 -o tara135_SRF_megahit_gather.csv \
    tara135_SRF_megahit.sig \
    mmetsp-k31-named.sbt.json \
    fungi/fungi-k31.sbt.json \
    protozoa/protozoa-k31.sbt.json \
    vertebrate_mammalian/vertebrate_mammalian-k31.sbt.json \
    vertebrate_other/vertebrate_other-k31.sbt.json \
    plant/plant-k31.sbt.json \
    invertebrate/invertebrate-k31.sbt.json
```

You should see an output that looks something like this:

```
overlap     p_query p_match avg_abund
---------   ------- ------- ---------
4.3 Mbp        3.3%    5.5%       1.0    Neoceratium fusus PA161109 MMETSP1074
1.2 Mbp        0.9%    3.8%       1.0    Calcidiscus leptoporus RCC1130 MMETSP...
4.2 Mbp        0.2%    0.3%       1.1    Neoceratium fusus PA161109 MMETSP1075
130.0 kbp      0.1%    0.2%       1.0    Brandtodinium nutricula RCC3387 MMETS...
90.0 kbp       0.1%    0.3%       1.0    Emiliania huxleyi CCMP370 MMETSP1157
70.0 kbp       0.0%    0.3%       1.0    Pelagomonas calceolata RCC969 MMETSP1328
70.0 kbp       0.1%    0.5%       1.2    Prasinoderma singularis RCC927 MMETSP...
80.0 kbp       0.0%    0.1%       1.2    Protoceratium reticulatum CCCM535=CCM...
50.0 kbp       0.0%    0.0%       1.2    Symbiodinium sp. CCMP421 MMETSP1110
50.0 kbp       0.0%    0.1%       1.0    Heterocapsa rotundata SCCAPK-0483 MME...
50.0 kbp       0.0%    0.0%       1.5    Scrippsiella trochoidea CCMP3099 MMET...
60.0 kbp       0.0%    0.1%       1.5    Prorocentrum minimum CCMP2233 MMETSP0269
found less than 40.0 kbp in common. => exiting

found 12 matches total;
the recovered matches hit 4.9% of the query
```

This will take about 5 to 10 minutes to run.

sourmash also outputs a csv with this and other information. We will use this csv
to visualize our results. 

We just ran sourmash on our assembly to give us an idea of the taxonomic 
make up of our sample in relation to known sequences. We could also run 
sourmash on our raw reads. This would give us an idea of the number of 
organisms we fail to assemble, or the amount of taxonomic diversity that
is missing from our assembly. Our sample has over 34 million reads and running
gather on a sample this size would take about an hour or two. We ran sourmash like
we did above, but with the reads, before this workshop. We will now download the
csv and compare the results from the assembly and from the reads. 

First, download the gather results from the raw reads. 

### TO DO:

+ finish interleaving reads.
+ run sourmash with scaled 10k on interleaved reads
+ post file to repo for download. 
+ write a script to compare the two

```
wget -O 
```


## Other notes

We will use the raw reads to do taxonomic classification. Although the errors
present in the raw reads may slow sourmash down a little, there is a chance 
when error trimming transcriptomes that real biological variation is removed. 
Because sourmash is so specific, the errors don't give false positives and so 
keeping all of the possible real variation could improve taxonomic recall in 
some cases.
