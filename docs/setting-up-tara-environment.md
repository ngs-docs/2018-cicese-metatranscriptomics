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
conda install fastqc multiqc trimmomatic khmer busco megahit sourmash salmon r dammit
```

To make sure you have access to these conda-installed programs, run the following:

```
echo "export PATH=/LUSTRE/apps/workshop/miniconda3/bin:$PATH" >> ~/.bashrc
```
and either 
  1. logout and log back in again, or
  2. source your bash profile to let the changes take effect: `source ~/.bashrc`

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

To put transrate in your path, you can execute:
```
echo "export PATH=/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64:$PATH" >> ~/.bashrc
source ~/.bashrc
```

Then check that transrate is properly installed with 
```
transrate -h
```


