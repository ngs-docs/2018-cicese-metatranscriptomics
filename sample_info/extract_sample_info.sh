

grep TARA_135 TARA_SAMPLES_CONTEXT_SEQUENCING_20170515.csv | grep Euk_DNA_RNA_ext | cut -d ',' -f 13,19,24 > TARA_135_samples.tsv
grep TARA_136 TARA_SAMPLES_CONTEXT_SEQUENCING_20170515.csv | grep Euk_DNA_RNA_ext | cut -d ',' -f 13,19,24 > TARA_136_samples.tsv
grep TARA_137 TARA_SAMPLES_CONTEXT_SEQUENCING_20170515.csv | grep Euk_DNA_RNA_ext | cut -d ',' -f 13,19,24 > TARA_137_samples.tsv

# in vim, with sed commands:
# %s/;$//g
# %s/,/<tab>/g
# %s/;/<tab>/g
#then use cut, paste (or awk) to make column #2 --> column #4 (sequencing type)

#if also want size fraction info, cut / paste from 1st column (sample name/info)

