# To edit these docs:

If you want to edit the docs, you'll need to grab the `mkdocs-material-dib` submodule. 

If cloning the repo for the first time:

```
git clone --recursive https://github.com/ngs-docs/2018-cicese-metatranscriptomics.git
```
`recursive` will pull the submodule as well as the main git repo.


If you already have the repo, and need to pull in the submodule:

```
git submodule update --init
```

## Editing and creating new docs

All docs should be written in md in the `docs` folder,
and committed to the repo as usual. If you added new 
documentation, be sure to add it into the site navigation 
in `mkdocs.yml`, or link to it from another doc.

 * sidebar navigation: edit `mkdocs.yml`
 
For the 2018 cicese workshop, also edit the schedule in `docs/index.md`


## After editing the docs, use mkdocs to deploy:

Install `mkdocs` and `ghp-import` if necessary:
```
conda install  -c conda-forge mkdocs
conda install -c conda-forge ghp-import
```

Build `docs` folder into `site`:

```
mkdocs build
```

View changes locally:
```
mkdocs serve
```

When satisfyed, push `docs` changes to `gh-pages` branch:

```
ghp-import site -p
```
