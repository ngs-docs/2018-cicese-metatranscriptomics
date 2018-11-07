# Evaluating a Metatranscriptome


## Generate assembly statistics with TransRate

We have installed transrate for you. However, if you need to install it in the future, the current installation instructions can be found [here](setting-up-tara-environment.md)

See options for running transrate:
```
/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64/transrate -h
```

Run transrate in `assembly` mode:
```
/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64/transrate --assembly
```









----


#### Assembly evaluation with rnaQUAST

rnaQUAST has some dependencies that conflict with other programs we have installed in the `tara` environment. To handle this, we'll create a new conda environment to run `rnaQUAST`

```
conda create -n rnaquast
source activate rnaquast
conda install rnaquast
```


```
python /LUSTRE/apps/workshop/miniconda3/envs/rnaquast/bin/rnaQUAST.py --transcripts final.contigs.fa -o rnaquast_results
```
