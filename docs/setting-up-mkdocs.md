# Setting up MkDocs in a new repo 

In your repo: 

Add `mkdocs-material-dib` as a submodule:

```
git submodule add https://github.com/dib-lab/mkdocs-material-dib.git
echo "site/" >> .gitignore
```

Commit:

```
git add .gitmodules .gitignore mkdocs-material-dib
git commit .gitmodules .gitignore mkdocs-material-dib -m 'Add mkdocs-material-dib submodule'
```

## Install `mkdocs` and `ghp-import` if necessary:

```
conda install  -c conda-forge mkdocs
conda install -c conda-forge ghp-import
```

## Edit site info and create docs

Edit the info in `mkdocs.yml` to reflect your repo name and address, color, and sidebar navigation. You'll also see that you can change the names for your docs directory and built web pages, but I assume here that they are `docs` and `site`, respectively.

Make a `docs` directory where you'll add the markdown docs for the website. Commit a simple `hello world` or similar md file to start. 

```
mkdir -p docs
echo '# Hello World' > docs/index.md
```

## Deploy the site

Go to `Settings` in your repo, and enable Github Pages for your repository.

In your repo, build the site:
```
mkdocs build
```

View your site locally:
```
mkdocs serve
```

Push your changes to `gh-pages`:
```
ghp-import site -p
```

Check your Github Pages URL for your simple md file:
```
https://<repo-owner>.github.io/<repo-name>
```

Continue to edit your docs, using `mkdocs build` to build the html version of the site, and `ghp-import` to push to `gh-pages`. Voila! 


## If you're cloning a repo with a mkdocs submodule:

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


For more details, follow the instructions [here](https://github.com/dib-lab/mkdocs-material-dib/tree/082e5399514cf2eb7c496eecb30a5570452966aa)
