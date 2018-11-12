# For Instructors:

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


## If setting up mkdocs for a new repo:

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

Edit the info in `mkdocs.yml` to reflect your repo name and address, color, and sidebar navigation. You'll also see that you can change the names for your docs directory and built web pages, but I assume here that they are `docs` and `site`, respectively.

Make a `docs` directory where you'll add the markdown docs for the website. Commit a simple `hello world` or similar md file to start. 

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

For more details, follow the instructions [here](https://github.com/dib-lab/mkdocs-material-dib/tree/082e5399514cf2eb7c496eecb30a5570452966aa)

