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
conda install fastqc multiqc trimmomatic khmer busco megahit sourmash salmon transrate-tools
#conda install dammit
```

Conda works for most of the software we'll use in the workshop. However, there are some exceptions>

To install transrate, follow the instructions [here](http://hibberdlab.com/transrate/installation.html).

```
cd <location to put transrate>
wget https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
tar xvf transrate-1.0.3-linux-x86_64.tar.gz
```

We installed transrate into the `/LUSTRE/apps/workshop/` directory, so it can be used as follows:

```
/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64/transrate -h 
```
