## What is conda?

Conda is a "package manager" or software installer. See [the full list of commands](https://conda.io/docs/commands.html).

`conda install` to install a package.

`conda list` to list installed packages.

`conda search` to search packages. Note that you'll see one package for every version of the software and for every version of Python (e.g. `conda search sourmash`).

## What is bioconda?

See [the bioconda paper](https://www.biorxiv.org/content/early/2017/10/27/207092) and the [bioconda web site](http://bioconda.github.io).

Bioconda is a community-enabled repository of 3,000+ bioinformatics packages, installable via the `conda` package
manager.  It consists of a set of recipes, [like this one, for sourmash](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/sourmash/meta.yaml), that are maintained   by the community.

## What problems does conda (and therefore bioconda) solve?

Conda tracks installed packages and their versions,
and makes sure that different installed packages don't have
conflicting dependencies. It makes software installation 
so much better!

### Constructing and using multiple environments

A feature that we do not use much here, but that can be very handy in some circumstances, is different environments.

"Environments" are multiple different collections of installed software. There are two reasons you might want to do this:

* first, you might want to try to exactly replicate a specific software install, so that you can replicate a paper or an old condition.
* second, you might be working with incompatible software, e.g. sometimes different software pipelines need different version of the same software. An example of this is older bioinformatics software that needs python2, while other software needs python3.


## Installing the software for this course

We already installed conda for you, but you will need to tell your terminal where
to look for that software. You can do this by executing the following:

```
echo 'export PATH=/LUSTRE/apps/workshop/miniconda3/bin:$PATH' >> ~/.bashrc
echo 'export PATH=/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64:$PATH' >> ~/.bashrc
source ~/.bashrc
```

After executing this, you should be able to run `conda` and look for software environments


To see the installed environments, run
```
conda info --envs
```

To activate the `tara` environment, which contains all software we'll use in the workshop, run

```
source activate tara
```

When you want to exit this environment later, you can execute `source deactivate` to return to the `base` env.


## Installing this software in the future

In the future, if you want to run the tutorials on your own, you'll
need to set up conda in your own account -
[see instructions](working-with-bioconda.md)).

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
echo 'export PATH=/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Then check that transrate is properly installed with 
```
transrate -h
```

## Finally --

Run 

```
~/works18
```

and then

```
source activate tara
```

so that we are all working on different distinct computers on the
`omica` cluster.
