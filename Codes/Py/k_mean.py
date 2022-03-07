# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:10 2021

@author: 12870
"""

import numpy as np
import pandas as pd
import os
import math
import re
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

os.chdir('D:\\project\\TCGA_stomach_analysis\\data')
rna = pd.read_csv('RNA_expr.csv', index_col=0)
val_name = rna.columns.tolist()
sample_l = []
for i in val_name:
    if re.match(r'^TCGA\.\w\w\.\w\w\w\w\.01A', i):
        i = 1
        sample_l.append(i)
    else:
        i = 0
        sample_l.append(i)
rna.columns = sample_l
rna = rna.T
rna.index.name = 'Type'
rna = rna.reset_index()

wgcna_gene = pd.read_csv('WGCNA_gene.csv')
name = wgcna_gene.iloc[:, 0].tolist()
gene_name = rna.columns.tolist()[1:]
for i in gene_name:
    if i in name:
        continue
    else:
        rna.drop(i, axis=1, inplace=True)
        
k = np.arange(1,11)
jarr = []
for i in k:
    model = KMeans(n_clusters=i)
    model.fit(X)
    jarr.append(model.inertia_)
    # 给这几个点打标
    plt.annotate(str(i),(i,model.inertia_))
plt.plot(k,jarr)
plt.show()

k = 4
model1 = KMeans(n_clusters=k)
# 跑模型
model1.fit(X)
# 需要知道每个类别有哪些参数
C_i = model1.predict(X)
# 还需要知道聚类中心的坐标
Muk = model1.cluster_centers_
plt.scatter(X[:,0],X[:,1],c=C_i,cmap=plt.cm.Paired)
# 画聚类中心
plt.scatter(Muk[:,0],Muk[:,1],marker='*',s=60)
for i in range(k):
    plt.annotate('center'+str(i + 1),(Muk[i,0],Muk[i,1]))
plt.show()
