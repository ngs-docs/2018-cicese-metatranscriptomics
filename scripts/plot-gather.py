#! /usr/bin/env python
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import pyupset as pyu
import matplotlib as mpl
from pickle import load
import pandas as pd
import glob 

genus_dict={}
for file in glob.glob('*csv'): 
    df=pd.read_csv(file, delimiter = ",")
    x=file.split('.')[0]
    genus_dict[x]=df
    print(x)

genus_dict['reads'] = genus_dict.pop('ERR1719497_paired_gather_all')
genus_dict['assembly'] = genus_dict.pop('tara_f135_full_megahit')

pplot=pyu.plot(genus_dict, unique_keys = ['name'])
pplot['figure'].savefig('plot-gather.png')
