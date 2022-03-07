# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 16:36:35 2021

@author: LanHao
"""
import re
import math


def dataProcessing(countsData, DEGnames):
    RNA_names = countsData.index.tolist()
    sample_name = countsData.columns.tolist()
    sample_name_l = []
    for gene in RNA_names:
        if gene in DEGnames:
            continue
        else:
            counts_data = countsData.drop(gene)
    for i in sample_name:
        if re.match(r'^TCGA\.(.*)\.(.*)\.01A\.(.*)\.(.*)', i):
            i = 1
            sample_name_l.append(i)
        else:
            i = 0
            sample_name_l.append(i)

    counts_data.columns = sample_name_l
    counts_data = counts_data.T
    counts_data.index.name = 'Type'
    counts_data.replace(0, 1, inplace=True)
    f = lambda num: math.log(num + 1, 10)
    counts_data = counts_data.applymap(f).reset_index()
    X = counts_data.iloc[:, 1:].values
    y = counts_data.iloc[:, 0].values
    return X, y, counts_data
