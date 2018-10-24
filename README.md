# 2018-metatranscriptomics-workshop-dev
(in progress) development repo for cicese metatranscriptomics workshop, nov 2018

To run the current pipeline on test data:

``` 
conda install snakemake
git clone https://github.com/bluegenes/2018-metatranscriptomics-workshop-dev
cd 2018-metatranscriptomics-workshop-dev
git clone https://github.com/bluegenes/rna_testdata.git # grab small rna test data
```

Test with a dry run:

```
snakemake --use-conda --configfile nemaRNA_config.yaml --dryrun
```

Actually run test data:

```
snakemake --use-conda --configfile nemaRNA_config.yaml
```

To run with tara data
 *note, you may need to change 'data_directory' within tara_config.yaml*
 *2nd note: ftp download is v. slow*

```
snakemake --use-conda --configfile tara_config.yaml
```




