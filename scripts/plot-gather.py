#! /usr/bin/env python
import pyupset as pyu
import matplotlib as mpl
import matplotlib.pyplot as plt
from pickle import load
import pandas as pd
import glob 

genus_dict={}
for file in glob.glob('*csv'): 
    df=pd.read_csv(file, delimiter = ",")
    x=file.split('.')[0]
    genus_dict[x]=df

genus_dict['reads'] = genus_dict.pop('ERR1719497_paired_gather_all')
genus_dict['assembly'] = genus_dict.pop('tara_f135_megahit_all_gather-scaled10k-k31')

pplot=pyu.plot(genus_dict, unique_keys = ['name'])
pplot['figure'].savefig('plot-gather.png')
