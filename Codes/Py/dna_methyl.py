# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 20:49:48 2021

@author: 12870
"""

import numpy as np
import pandas as pd
import os
import re
import math
from collections import Counter

os.chdir('D:\\project\\TCGA_stomach_analysis\\data')

met = pd.read_csv('HumanMethylation450.tsv', sep='\t', index_col=0)
met.dropna(axis=0, how='any', inplace=True)
met = met.sort_index()
map_df = pd.read_csv('illuminaMethyl450_hg38_GDC.tsv',
                     sep='\t',
                     ).iloc[:, 0:2]
map_df.sort_index(inplace=True)
map_df.replace('.', np.NAN, inplace=True)
map_df.dropna(axis=0, how='any', inplace=True)
cpg_list = map_df.index.tolist()
met.reset_index(inplace=True)
l = met.iloc[:, 0].tolist()
met['id']=l
met.set_index('id', inplace=True)
for i in met.index.tolist():
    if i in cpg_list:
        met.loc[i, 'sample'] = map_df.loc[i, 'gene']


map_df.set_index('gene')