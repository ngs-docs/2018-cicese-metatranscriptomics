# 2018-metatranscriptomics-workshop-dev
(in progress) development repo for cicese metatranscriptomics workshop, nov 2018

to run the current pipeline on test data:

conda install snakemake

git clone https://github.com/bluegenes/2018-metatranscriptomics-workshop-dev
cd 2018-metatranscriptomics-workshop-dev
git clone https://github.com/bluegenes/rna_testdata.git # grab small rna test data

# test with a dry run:
snakemake --use-conda --configfile nemaRNA_config.yaml --dryrun

#actually run:
snakemake --use-conda --configfile nemaRNA_config.yaml



