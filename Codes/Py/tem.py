# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:06:22 2021

@author: 12870
"""

import numpy as np
import pandas as pd
import os
import re
import math
from collections import Counter

os.chdir('D:\\project\\TCGA_stomach_analysis\\data')
a = pd.read_csv('./download_rnaseq/mRNA/mRNA_expr.csv',
                index_col=0)
b = pd.read_csv('./download_mirna/miRNA_expr.csv',
                index_col=0)
a_ = a.columns.tolist()
b_ = b.columns.tolist()
for i in range(len(a_)):
    a_[i] = re.sub(r'\.[0-9][0-9]R\S+', '', a_[i])
for i in range(len(b_)):
    b_[i] = re.sub(r'\.[0-9][0-9]R\S+', '', b_[i])
val_name = list(set(a_) & set(b_))
c = list(set(a_) ^ set(b_))

sample_name_l = []
for i in val_name:
    if re.match(r'^TCGA\.\w\w\.\w\w\w\w\.01A', i):
        i = 1
        sample_name_l.append(i)
    else:
        i = 0
        sample_name_l.append(i)
Counter(sample_name_l)

sample_name = []
for i in c:
    if re.match(r'^TCGA\.\w\w\.\w\w\w\w\.01A', i):
        i = 1
        sample_name.append(i)
    else:
        i = 0
        sample_name.append(i)
Counter(sample_name)

mRNA = pd.read_csv('./download_rnaseq/mRNA/mRNA_expr.csv',
                   index_col=0)
mrna_sample_name = mRNA.columns.tolist()
for i in range(len(mrna_sample_name)):
    mrna_sample_name[i] = re.sub(r'\.[0-9][0-9]R\S+',
                                 '',
                                 mrna_sample_name[i])
mRNA.columns = mrna_sample_name
for i in mrna_sample_name:
    if i in val_name:
        continue
    else:
        mRNA.drop(i, axis=1, inplace=True)
mRNA_DEG_name = list(pd.read_csv('mRNA_DEGname.csv').iloc[:, 0])
mrna_gene_name = mRNA.index.tolist()
for i in range(len(mrna_gene_name)):
    mrna_gene_name[i] = 'gene-' + mrna_gene_name[i]
mRNA.index = mrna_gene_name
mRNA = mRNA.T
for i in mrna_gene_name:
    if i in mRNA_DEG_name:
        continue
    else:
        mRNA.drop(i, axis=1, inplace=True)


lncRNA = pd.read_csv('./download_rnaseq/lncRNA/lncRNA_expr.csv',
                     index_col=0)
lncrna_sample_name = lncRNA.columns.tolist()
for i in range(len(lncrna_sample_name)):
    lncrna_sample_name[i] = re.sub(r'\.[0-9][0-9]R\S+',
                                   '',
                                   lncrna_sample_name[i])
lncRNA.columns = lncrna_sample_name
for i in lncrna_sample_name:
    if i in val_name:
        continue
    else:
        lncRNA.drop(i, axis=1, inplace=True)
lncRNA_DEG_name = list(pd.read_csv('lncRNA_DEGname.csv').iloc[:, 0])
lncrna_gene_name = lncRNA.index.tolist()
for i in range(len(lncrna_gene_name)):
    lncrna_gene_name[i] = 'lnc-' + lncrna_gene_name[i]
lncRNA.index = lncrna_gene_name
lncRNA = lncRNA.T
for i in lncrna_gene_name:
    if i in lncRNA_DEG_name:
        continue
    else:
        lncRNA.drop(i, axis=1, inplace=True)


miRNA = pd.read_csv('./download_mirna/miRNA_expr.csv',
                    index_col=0)
mirna_sample_name = miRNA.columns.tolist()
for i in range(len(mirna_sample_name)):
    mirna_sample_name[i] = re.sub(r'\.[0-9][0-9]R\S+',
                                  '',
                                  mirna_sample_name[i])
miRNA.to_csv('miRNA.csv')
miRNA = pd.read_csv('miRNA.csv', index_col=0)
miRNA.columns = mirna_sample_name
for i in list(set(mirna_sample_name)):
    if i in val_name:
        continue
    else:
        miRNA.drop(i, axis=1, inplace=True)
miRNA_DEG_name = list(pd.read_csv('miRNA_DEGname.csv').iloc[:, 0])
mirna_gene_name = miRNA.index.tolist()
miRNA = miRNA.T
for i in mirna_gene_name:
    if i in miRNA_DEG_name:
        continue
    else:
        miRNA.drop(i, axis=1, inplace=True)
RNA_expr = pd.concat((mRNA, lncRNA, miRNA), axis=1)
RNA_expr.to_csv('RNA_expr.csv')
RNA_expr.index.name = 'Type'
RNA_expr.replace(0, 1, inplace=True)
f = lambda num: math.log(num + 1, 10)
RNA_expr = RNA_expr.applymap(f)
RNA_expr = RNA_expr.T
sample_attributes = pd.DataFrame(index=val_name)
sample_l = []
for i in val_name:
    if re.match(r'^TCGA\.\w\w\.\w\w\w\w\.01A', i):
        i = 0
        sample_l.append(i)
    else:
        i = 1
        sample_l.append(i)
sample_attributes['Normal'] = sample_l
sample_attributes.to_csv('attributes.csv')


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







clin = pd.read_csv('clinical1.csv')
clin_sample_name = clin['submitter_id'].tolist()
for i in range(len(clin_sample_name)):
    clin_sample_name[i] = re.sub('-', '.', clin_sample_name[i])
RNA = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\rna_expr.csv',
                  index_col=0)
rna_sample_name = RNA.columns.tolist()
for i in range(len(rna_sample_name)):
    rna_sample_name[i] = re.sub(r'\.[0-1][0-1][A-Z]$', '', rna_sample_name[i])
RNA.columns = rna_sample_name
common_sample_name = []
common_sample_name = list(set(rna_sample_name) & set(clin_sample_name))
for i in rna_sample_name:
    if i in clin_sample_name:
        common_sample_name.append(i)
for i in rna_sample_name:
    if i in common_sample_name:
        continue
    elif i in RNA.columns.tolist():
        RNA.drop(i, axis=1, inplace=True)
a = pd.read_csv('after_rfe.csv')
rna_n = a.columns.tolist()[1:]
rna_n_ = RNA.index.tolist()
for i in rna_n_:
    if i in rna_n:
        continue
    else:
        RNA.drop(i, axis=0, inplace=True)
RNA = RNA.T
rna.to_csv('rna_clinical.csv')

rna = pd.read_csv('rna_clinical.csv', index_col=0)
rna = rna.T
clinical = pd.read_excel('clinical1.xlsx', index_col=0)
clinical.drop(clinical[clinical.vital_status == 'Not Reported'].index,
              axis=0,
              inplace=True)
clinical.dropna(axis=0,
                how='all',
                subset=['days_to_last_follow_up', 'days_to_death'],
                inplace=True)
clinical = clinical.replace('nan', np.NAN)
days_to_death = clinical['days_to_death'].tolist()
days_to_last_follow_up = clinical['days_to_last_follow_up'].tolist()
for i in range(clinical.shape[0]):
    if pd.isna(pd.Series(days_to_last_follow_up)[i]):
        days_to_last_follow_up[i] = days_to_death[i]
clinical['days_to_last_follow_up'] = days_to_last_follow_up
clinical.rename(columns={'days_to_last_follow_up': 'OS_time'}, inplace=True)
clinical.drop('days_to_death', axis=1, inplace=True)
clinical.index.name = 'id'
clinical_name = clinical.index.tolist()
sample_name = rna.index.tolist()
for i in clinical_name:
    if i in sample_name:
        continue
    elif i in clinical.index.tolist():
        clinical.drop(i, axis=1, inplace=True)
a = pd.DataFrame([dict(Counter(rna.columns.tolist()))]).T
a.drop(a[a.iloc[:, 0] == 1].index, inplace=True)
b = a.index.tolist()
rna = rna.T
for i in b:
    rna[i].mean(axis=1)
new_name = [i + '-mean' for i in b ]
for i in new_name:
    for j in b:
        if j + '-mean' == i:
            rna[i] = rna[[j,j]].mean(axis=1)
for i in b:
    rna.drop(i, axis=1, inplace=True)
for i in range(len(new_name)):
    new_name[i] = new_name[i].split('-')[0]
sample_n = rna.columns.tolist()
sample_n[-26:] = new_name
rna.columns = sample_n
rna = rna.T
rna.index.name = 'id'
common = list(set(sample_n) & set(clinical_name))
for i in sample_n:
    if i in common:
        continue
    else:
        rna.drop(i, axis=0, inplace=True)
for i in clinical_name:
    if i in common:
        continue
    else:
        clinical.drop(i, axis=0, inplace=True)
rna = rna.sort_index()
clinical = clinical.sort_index()
cox_km = pd.concat((rna, clinical), axis=1)
cox_km.to_csv('cox_km.csv')
cox_km = pd.read_csv('cox_km.csv', index_col=0)
cox_km['risk_score'] = np.zeros(cox_km.shape[0])
for i in range(cox_km.shape[0]):
    cox_km.iloc[i, -1] = np.sum(np.array(cox_km.iloc[i, :-3]) * np.array(cox_km.iloc[349, :-3]))
cox_km.to_csv("multi_cox_3mRNA_risk_score_DFS.csv")
cox_km.drop('coef', axis=0, inplace=True)
risk = list(np.array(cox_km.iloc[:, -1]))
med = np.median(risk)
for i in range(len(risk)):
    if i == 0:
        risk[i] = 'low risk'
    else:
        risk[i] = 'high risk'
cox_km['Risk'] = risk
cox_km.drop(cox_km[cox_km.OS_time == 0].index, axis=0, inplace=True)

rna = pd.read_csv('rna_expr.csv', index_col=0)
after = pd.read_csv('after_rfe.csv', index_col=0)
gene = after.columns.tolist()
rna_gene = rna.index.tolist()
for i in rna_gene:
    if i in gene:
        continue
    else:
        rna.drop(i, axis=0, inplace=True)
rna.to_csv('tme_rna.csv')

a = pd.read_csv('expr_cpm.csv', index_col=0)
name = a.index.tolist()
for i in range(0, 2060):
    name[i] = name[i].split('-', 1)
    name[i] = name[i][1]
a.index=name
a.to_csv('tme_rna.csv')


b = pd.read_excel('score.xlsx')
c = b.iloc[:, 0].tolist()
l = []
for i in c:
    if re.match(r'^TCGA\.(.*)\.(.*)\.01A\.(.*)\.(.*)', i):
        
        l.append('Tumor')
    else:
        
        l.append('Normal')
b['Type'] = l
b.to_csv('score_t.csv')
