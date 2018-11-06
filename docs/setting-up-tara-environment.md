# Installing all software for the workshop

After setting up conda ([here](working-with-bioconda.md)), let's install some software.

First, let's create an environment we'll use for the course:

``` 
conda create -n tara
```

and activate it:

```
source activate tara
```

Now, let's install the software!

```
conda install fastqc multiqc trimmomatic khmer busco megahit sourmash salmon
```
