## Installing this software in the future

In the future, if you want to run the tutorials on your own, you'll
need to set up conda in your own account -
[see instructions](working-with-bioconda.md)).

To install all the tools we used in the future, we can 1. install from
an environment `yaml` file, or 2. install everything manually.

### Option 1: Installing from a file

At the end of the course, we exported info from the `tara` environment using:
```
conda env export -n tara -f $PROJECT/tara_conda_environment.yaml
```
You can find download this file [here](files/tara_conda_environment.yaml)


To make a new conda env like this, download that file, then run:

```
conda env create -f tara_conda_environment.yaml
```

you shuld then be able to:

```
source activate tara
```

and you remember you can exit this environment with `source deactivate`


### Option 2: Installing Manualling


Then create an environment to work in:

``` 
conda create -n tara
```

and activate it:

```
source activate tara
```

Now, install the software!

```
conda install fastqc multiqc trimmomatic khmer busco megahit sourmash salmon r dammit cd-hit
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
echo 'export PATH=/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Then check that transrate is properly installed with 
```
transrate -h
```

