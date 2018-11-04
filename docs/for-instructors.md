# For Instructors:

To use mkdocs to deploy lessons:

#### If mkdocs is already set up

All docs should be written in md in the `docs` folder,
and committed to the repo as usual.

Install `mkdocs` and `ghp-import` if necessary:
```
conda install  -c conda-forge mkdocs
conda install -c conda-forge ghp-import
```

Build docs folder into site/
```
mkdocs build to update the docs
```

View changes locally:
```
mkdocs serve
```

Push docs changes to gh-pages branch
```
ghp-import site -p
```


#### Inital Setup 

Grab the mkdocs repo, dib lab flavor, and follow setup instructions
```
git clone https://github.com/dib-lab/mkdocs-material-dib.git
```

