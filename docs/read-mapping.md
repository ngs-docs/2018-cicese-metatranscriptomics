# Mapping reads to a metatranscriptome




## Quantification with Salmon

_testing_

## Index the Assembly

```
#for assembly in *_megahit.fasta

for assembly in *.contigs.fa
do
  base=$(basename $assembly .contigs.fa)
  echo $base

  salmon index -t ${base}.contigs.fa -i ${base}_salmon --threads 2
  #salmon index -t ${base}.fasta -i ${base}_salmon --threads 2
done
```

## Run quantification

```
for filename in *_1.qc.fq.gz
do
  base=$(basename $filename _1.qc.fq.gz)
  echo $base
  
  for assembly_index in *_salmon
    do
    time salmon quant -i ${assembly_index} -l A -1 ${base}_1.qc.fq.gz -2 ${base}_2.qc.fq.gz -o ${base}_quant  -p 2
    done
done
```


