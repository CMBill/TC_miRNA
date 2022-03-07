# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 16:40:45 2021

@author: LanHao
"""

import pandas as pd


def readFiles(filepath_counts, filepath_names):
    counts_data = pd.read_csv(filepath_counts, sep=',', index_col=0)
    DEGnames = list(pd.read_csv(filepath_names,
                                sep=',',
                                index_col=False).iloc[:, 0])
    return counts_data, DEGnames
