# Software Installation with Bioconda 

Learning objectives:  

  * learn what bioconda is
  * understand basic `conda` commands
  * learn how to list installed software packages 
  * learn how to manage multiple installation environments

## What is bioconda?

See [the bioconda paper](https://www.biorxiv.org/content/early/2017/10/27/207092) and the [bioconda web site](http://bioconda.github.io).

Bioconda is a community-enabled repository of 3,000+ bioinformatics packages, installable via the `conda` package
manager.  It consists of a set of recipes, [like this one, for sourmash](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/sourmash/meta.yaml), that are maintained by the community.

It just works, and it's effin' magic!!

## What problems does conda (and therefore bioconda) solve?

Conda tracks installed packages and their versions, 
and makes sure that different installed packages don't have
conflicting dependencies (we'll explain this below).

## Installing conda and enabling bioconda

Download and install miniconda in your home directory:
```
cd 
curl -O -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3 #install in $HOME directory
echo export PATH=$PATH:/$HOME/miniconda3/bin >> ~/.bashrc
```

Then, run the following command (or start a new terminal session) in order to activate the conda environment:

```
source ~/.bashrc
```

Configure channels for installing software.
Note: It is important to add them in this order so that 
the priority is set correctly (that is, conda-forge is highest priority).
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Try installing something:

```
conda install sourmash
```

and running it --
```
sourmash
```
will produce some output. (We'll tell you more about sourmash later.)

yay!

## Using conda

Conda is a "package manager" or software installer. See [the full list of commands](https://conda.io/docs/commands.html).

`conda install` to install a package.

`conda list` to list installed packages.

`conda search` to search packages. Note that you'll see one package for every version of the software and for every version of Python (e.g. `conda search sourmash`).

During this workshop, we have been using the `tara` environment. We can expor this environment into a "yaml" text file to keep track 
of everything we have installed, including versions!

```
`conda env export -n tara -f $PROJECT/tara_conda_environment.yaml`
```

You can read more about exporting environments [here](https://conda.io/docs/commands/env/conda-env-export.html)

## Using bioconda

bioconda is a channel for conda, which just means that you
can "add" it to conda as a source of packages. That's what the `conda config` above does.

Note, Bioconda supports only 64-bit Linux and Mac OSX.

You can check out [the bioconda site](https://bioconda.github.io/).

### Finding bioconda packages

You can use `conda search`, or you can use google, or you can go visit [the list of recipes](https://bioconda.github.io/recipes.html#recipes).

### Constructing and using multiple environments

A feature that we do not use much here, but that can be very handy in some circumstances, is different environments.

"Environments" are multiple different collections of installed software. There are two reasons you might want to do this:

* first, you might want to try to exactly replicate a specific software install, so that you can replicate a paper or an old condition.
* second, you might be working with incompatible software, e.g. sometimes different software pipelines need different version of the same software. An example of this is older bioinformatics software that needs python2, while other software needs python3.

To create a new environment named `pony`, type:

```
conda create -n pony
```

Then to activate (switch to) that environment, type:

```
source activate pony
```

And now when you run `conda install`, it will install packages into this new environment, e.g.
```
conda install -y checkm-genome
```
(note here that checkm-genome *requires* python 2).


### Freezing an environment

This will save the list of **conda-installed** software you have in a particular
environment to the file `packages.txt`:

```
conda list --export packages.txt
```
(it will not record the software versions for software not installed by conda.)

```
conda install --file=packages.txt
```
will install those packages in your local environment.
To list environments, type:
```
conda env list
```
and you will see that you have two environments, `base` and
`pony`, and `pony` has a `*` next to it because that's your
current environment.

### Leaving an environment

And finally, to switch back to your base environment, do:

```
source deactivate
```
and you'll be back in the original environment.

### Meditations on reproducibility and provenance

If you want to impress reviewers and also keep track of
what your software versions are, you can:

* manage all your software inside of conda
* use `conda list --export software.txt` to create a list of all your software and put it in your supplementary material.

This is also something that you can record for yourself, so that if you are trying to exactly reproduce 

### Using it on your own compute system (laptop or HPC)

conda works on Windows, Mac, and Linux.

bioconda works on Mac and Linux.

It does not require admin privileges to install, so you can
install it on your own local cluster quite easily.
