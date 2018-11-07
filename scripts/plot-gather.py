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

pplot=pyu.plot(genus_dict, unique_keys = ['name'])
pplot['figure'].savefig('test.png')
