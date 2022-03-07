# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 16:41:52 2021

@author: 12870
"""
import pandas as pd
import numpy as np
import math
import re
data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\expr.csv',
                   sep=',',
                   index_col=0)
name = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\gene_name_wgcna.xlsx',
                     index_col=False)
RNA_names = data.index.tolist()
name = name.iloc[:, 0].tolist()
sample_name = data.columns.tolist()
sample_name_l = []

for gene in RNA_names:
    if gene in name:
        continue
    else:
        data.drop(gene, axis=0, inplace=True)
for i in sample_name:
    if re.match(r'^TCGA\.(.*)\.(.*)\.01A\.(.*)\.(.*)', i):
        i = 1
        sample_name_l.append(i)
    else:
        i = 0
        sample_name_l.append(i)

data.columns = sample_name_l
data = data.T
data.index.name = 'Type'
data.replace(0, 1, inplace=True)
f = lambda num: math.log(num + 1, 10)
data = data.applymap(f).reset_index()
data.to_csv('D:\\project\\TCGA_stomach_analysis\\data\\expr_after_first_selection.csv')
a = data.sum(axis=0)
for i in sample_name:
    data[i] = data[i] / a[i] * 1000000
    
f = lambda num: math.log(num + 1, 10)
data = data.applymap(f)
data.to_csv('D:\\project\\TCGA_stomach_analysis\\data\\expr_cpm.csv')

clin_info = data.columns.tolist()
cliInfo = pd.DataFrame(index=clin_info, columns=['Normal', 'Tumor'])

for i in clin_info:
    if re.match(r'^TCGA\.(.*)\.(.*)\.01A\.(.*)\.(.*)', i):
        cliInfo.loc[i, 'Normal'] = 0
        cliInfo.loc[i, 'Tumor'] = 1
    else:
        cliInfo.loc[i, 'Normal'] = 1
        cliInfo.loc[i, 'Tumor'] = 0
cliInfo.to_csv('D:\\project\\TCGA_stomach_analysis\\data\\cliInfo.csv')





data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\expr.csv',
                   sep=',',
                   index_col=0)
X, y = np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0])
from collections import Counter
from sklearn.datasets import make_classification
X, y = make_classification(n_classes=2, class_sep=2,
                           weights=[0.9, 0.1], n_informative=3, 
                           n_redundant=1, flip_y=0,
                           n_features=20, n_clusters_per_class=1, 
                           n_samples=1000, random_state=10)

from imblearn.over_sampling import SMOTE
smo = SMOTE(random_state=42)
X_smo, y_smo = smo.fit_resample(X, y)
print(Counter(y_smo))

data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\TOM.csv',
                   sep=',')
name = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\RNA_DEGname.xlsx',
                   index_col=False)
name = name.iloc[:, 0].tolist()
data.columns=name
data.index = name
name_905 = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\gene_name_wgcna.xlsx',
                     index_col=False)
name_905 = name_905.iloc[:, 0].tolist()
data.loc['gene-ENPP7','gene-ENPP7']
for i in name:
    if i in name_905:
        continue
    else:
        data.drop(i, axis=0, inplace=True)
        data.drop(i, axis=1, inplace=True)
data.to_csv('D:\\project\\TCGA_stomach_analysis\\data\\TOM.csv')

data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\expr_after_first_selection.csv')
sample = data.columns.tolist()
for i in sample:
    if i in name:
name = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\GENERANK100.xlsx')
name = name.iloc[:, 0].tolist()
        continue
    else:
        data.drop(i, axis=1, inplace=True )
data.to_csv('D:\\project\\TCGA_stomach_analysis\\data\\generank_selected.csv')




mrna = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\download_rnaseq\\mRNA\\mRNA_expr.csv',
                   index_col=0)
sample_name_m = mrna.columns.tolist()
mirna = pd.read_csv('D:\project\TCGA_stomach_analysis\data\download_mirna\miRNA_expr1.csv',
                    index_col=0)
sample_name_mi = mirna.columns.tolist()
for i in range(len(sample_name_mi)):
    sample_name_mi[i] = re.sub(r'.13$', '', sample_name_mi[i])
mirna.columns =sample_name_mi
val_name = list(set(sample_name_m)^set(test_name))
test_name = list(set(sample_name_m)&set(sample_name_mi))
name = list(set(sample_name_mi) - set(test_name))
for i in sample_name_m:
    if i in val_name:
        continue
    else:
        mrna.drop(i, axis=1, inplace=True)
        
lncrna=pd.read_csv('D:\project\TCGA_stomach_analysis\data\download_rnaseq\lncRNA\lncRNA_expr.csv',
                   index_col=0)
sample_name_lnc = lncrna.columns.tolist()
for i in sample_name_lnc:
    if i in val_name:
        continue
    else:
        lncrna.drop(i, axis=1, inplace=True)
lnc_gene=lncrna.index.tolist()
for i in range(len(lnc_gene)):
    lnc_gene[i] = 'lnc-'+lnc_gene[i]
lncrna.index=lnc_gene
m_gene = mrna.index.tolist()
for i in range(len(m_gene)):
    m_gene[i] = 'gene-'+m_gene[i]
mrna.index=m_gene

for i in sample_name_mi:
    if i in name:
        continue
    else:
        mirna.drop(i, axis=1, inplace=True)
rna_val = pd.concat((mrna, lncrna), axis=0)
sample_name_l=[]
for i in val_name:
    if re.match(r'^TCGA\.(.*)\.(.*)\.01A\.(.*)\.(.*)', i):
        i = 1
        sample_name_l.append(i)
    else:
        i = 0
        sample_name_l.append(i)
rna_val.columns = sample_name_l
rna_val = rna_val.T
rna_val.index.name = 'Type'
rna_val.replace(0, 1, inplace=True)
f = lambda num: math.log(num + 1, 10)
rna_val = rna_val.applymap(f).reset_index()
rna_val_name = rna_val.columns.tolist()
del(rna_val_name[0])
for i in rna_val_name:
    if i in ganame:
        continue
    else:
        rna_val.drop(i, axis=1, inplace=True)

ganame = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\GAname.xlsx')
ganame = ganame.iloc[:, 0].tolist()
lncrna_name = lncrna.columns.tolist()
for i in range(1,len(lncrna_name)):
    lncrna_name[i] = 'lnc-'+lncrna_name[i]
lncrna.columns = lncrna_name
rna_val.to_csv('D:\\project\\TCGA_stomach_analysis\\data\\valirna.csv')

mrna.columns = sample_name_l
mrna = mrna.T
mrna.index.name = 'Type'
mrna.replace(0, 1, inplace=True)
f = lambda num: math.log(num + 1, 10)
mrna = mrna.applymap(f).reset_index()
ganame = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\GAname.xlsx')
ganame = ganame.iloc[:, 0].tolist()
mrna_name = mrna.columns.tolist()
for i in range(1,len(mrna_name)):
    mrna_name[i] = 'gene-'+mrna_name[i]
mrna.columns = mrna_name
del(mrna_name[0])
for i in mrna_name:
    if i in ganame:
        continue
    else:
        mrna.drop(i, axis=1, inplace=True)
mrna.to_csv('D:\\project\\TCGA_stomach_analysis\\data\\valimrna.csv')

n = data.columns.tolist()
m = rna_val.columns.tolist()
del(n[0])
del(m[0])
for i in m:
    if i in n:
        continue
    else:
        rna_val.drop(i, axis=1, inplace=True)
rna_val.to_csv('iitest_9_32.csv')
