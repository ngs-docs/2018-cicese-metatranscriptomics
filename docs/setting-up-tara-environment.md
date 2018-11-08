# Installing all software for the workshop

We already installed all course software for you, but you need to tell your terminal where
to look for that software. You can do this by executing the following:

```
echo "export PATH=/LUSTRE/apps/workshop/miniconda3/bin:$PATH" >> ~/.bashrc
echo "export PATH=/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64:$PATH" >> ~/.bashrc
source ~/.bashrc
```


# Installing software in the future

First, make sure you've set up conda ([here](working-with-bioconda.md)).

Then  create an environment to work in:

``` 
conda create -n tara
```

and activate it:

```
source activate tara
```

Now, install the software!

```
conda install fastqc multiqc trimmomatic khmer busco megahit sourmash salmon r dammit
```

To make sure you have access to these conda-installed programs, try to execute some of them:

```
sourmash
```

Conda works for most of the software we'll use in the workshop. However, there are some exceptions: notably, transrate. 

To install transrate, follow the instructions [here](http://hibberdlab.com/transrate/installation.html).

You can also try the following:
```
cd <location-to-put-transrate>
wget https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
tar xvf transrate-1.0.3-linux-x86_64.tar.gz
```

To put transrate in your path, you can execute:
```
echo "export PATH=<path-to-transrate>/transrate-1.0.3-linux-x86_64:$PATH" >> ~/.bashrc
source ~/.bashrc
```

Then check that transrate is properly installed with 
```
transrate -h
```


