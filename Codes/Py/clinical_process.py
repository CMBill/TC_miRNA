# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 15:12:30 2021

@author: 12870
"""

import numpy as np
import pandas as pd
import os
import re
import math
from collections import Counter

os.chdir('D:\\project\\TCGA_stomach_analysis\\data')

# 准备好临床数据，准备合并
clin = pd.read_csv('clinical.csv')
Counter(list(clin['ajcc_pathologic_stage']))
Counter(list(clin['vital_status']))
d = dict(Counter(list(clin['days_to_last_follow_up'])))

clin.dropna(axis=0,
            how='any',
            subset=['ajcc_pathologic_stage','age_at_index'],
            inplace=True)
clin.drop(clin[(clin.iloc[:, 2] < 30) & (clin.iloc[:, 2] >= 0)].index,
          axis=0,
          inplace=True)
clin = clin.replace('nan', np.NAN)
days_to_death = clin['days_to_death'].tolist()
days_to_last_follow_up = clin['days_to_last_follow_up'].tolist()
for i in range(clin.shape[0]):
    if pd.isna(pd.Series(days_to_last_follow_up)[i]):
        days_to_last_follow_up[i] = days_to_death[i]
clin['days_to_last_follow_up'] = days_to_last_follow_up
clin.rename(columns={'days_to_last_follow_up': 'OS_time'}, inplace=True)
clin.drop('days_to_death', axis=1, inplace=True)

clin.drop([344, 349, 387], axis=0, inplace=True)
clin.drop(clin[clin.vital_status == 'Not Reported'].index,
              axis=0,
              inplace=True)
clin['vital_status'].replace('Alive', 0, inplace=True)
clin['vital_status'].replace('Dead', 1, inplace=True)
clin['gender'].replace('female', 0, inplace=True)
clin['gender'].replace('male', 1, inplace=True)
clin['ajcc_pathologic_stage'].replace(['Stage I', 'Stage IB', 'Stage IA'],
                                      0,
                                      inplace=True)
clin['ajcc_pathologic_stage'].replace(['Stage II', 'Stage IIA', 'Stage IIB'],
                                      1,
                                      inplace=True)
clin['ajcc_pathologic_stage'].replace(['Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC'],
                                      2,
                                      inplace=True)
clin['ajcc_pathologic_stage'].replace(['Stage IV'],
                                      3,
                                      inplace=True)
clin.drop(clin[(clin.ajcc_pathologic_t == 'TX')].index,
          axis=0,
          inplace=True)
clin.drop(clin[(clin.ajcc_pathologic_m == 'MX')].index,
          axis=0,
          inplace=True)
clin.drop(clin[(clin.ajcc_pathologic_n == 'NX')].index,
          axis=0,
          inplace=True)
clin['ajcc_pathologic_t'].replace(['T1', 'T1a', 'T1b'],
                                  0,
                                  inplace=True)
clin['ajcc_pathologic_t'].replace(['T2', 'T2a', 'T2b'],
                                  1,
                                  inplace=True)
clin['ajcc_pathologic_t'].replace('T3',
                                  1,
                                  inplace=True)
clin['ajcc_pathologic_t'].replace(['T4', 'T4a', 'T4b'],
                                  1,
                                  inplace=True)

clin['ajcc_pathologic_n'].replace('N0', 
                                  0,
                                  inplace=True)
clin['ajcc_pathologic_n'].replace('N1',
                                  1,
                                  inplace=True)
clin['ajcc_pathologic_n'].replace('N2',
                                  2,
                                  inplace=True)
clin['ajcc_pathologic_n'].replace(['N3', 'N3a', 'N3b'],
                                  3,
                                  inplace=True)
clin['ajcc_pathologic_m'].replace('M0',
                                  0,
                                  inplace=True)
clin['ajcc_pathologic_m'].replace('M1',
                                  1,
                                  inplace=True)
clin.dropna(axis=0,
            how='all',
            subset=['ajcc_pathologic_n'],
            inplace=True)
clin = clin.set_index('submitter_id')
clin.index.name = 'id'
clin['OS_time'] = np.array(clin['OS_time']) / 30
age = np.array(clin['age_at_index'])
age = np.where(age < 50, 0, age)
age = np.where((age < 60) & (age >= 50), 1, age)
age = np.where((age < 70) & (age >= 60), 2, age)
age = np.where((age < 80) & (age >= 70), 3, age)
age = np.where((age <= 90) & (age >= 80), 4, age)
clin['age_at_index'] = age
clin_sample_name = clin.index.tolist()

# 准备挑选出的基因的数据
rna = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\rna_expr.csv',
                  index_col=0)

a = pd.read_csv('after_rfe.csv')
rna_n = a.columns.tolist()[1:]
rna_n_ = rna.index.tolist()
for i in rna_n_:
    if i in rna_n:
        continue
    else:
        rna.drop(i, axis=0, inplace=True)
rna_sample_name = rna.columns.tolist()

# 删除正常样本
for i in rna_sample_name:
    if re.match(r'^TCGA\.\w\w\.\w\w\w\w\.01A', i):
        continue
    else:
        rna.drop(i, axis=1, inplace=True)
rna_sample_name = rna.columns.tolist()
for i in range(len(rna_sample_name)):
    rna_sample_name[i] = re.sub(r'\.[0-1][0-1][A-Z]$', '', rna_sample_name[i])
    rna_sample_name[i] = re.sub('\.', '-', rna_sample_name[i])
rna.columns = rna_sample_name

#挑选共同样本，合并数据
len(set(rna_sample_name))
len(set(clin_sample_name))
common_name = list(set(rna_sample_name) & set(clin_sample_name))

for i in rna_sample_name:
    if i in common_name:
        continue
    else:
        rna.drop(i, axis=1, inplace=True)
for j in clin_sample_name:
    if j in common_name:
        continue
    else:
        clin.drop(j, axis=0, inplace=True)

rna = rna.T
clin = clin.sort_index()
rna = rna.sort_index()
cox_km = pd.concat((rna, clin), axis=1)
cox_km.to_csv('cox_km_c.csv')

cox_km_c = pd.read_csv('cox_km_c.csv', index_col=0)
l_gene = ['gene-ADAMTS18', 'gene-AKR1B15', 'gene-APOC1', 'gene-CAMKV',
          'gene-PPBP', 'gene-CHRNA1', 'gene-COL9A3', 'gene-DNASE1L3',
          'gene-HOXC8', 'lnc-KCNMB2-AS1']
gene_uni = cox_km.columns.tolist()[0:39]
for i in gene_uni:
    if i in l_gene:
        continue
    else:
        cox_km_c.drop(i, axis=1, inplace=True)
cox_km_c.to_csv('multi_1.csv')

# 风险评分数据准备
cox_km_r = pd.read_csv('multi_1.csv')
coef = np.array([0.7111, 0.4144, 0.6337, 0.4154,
                 0.3643, -0.5461, -0.4738, -0.5615,
                 -0.4286, 0.4436])
gene_ex =np.array(cox_km_r.iloc[:, 1:(len(coef) + 1)])
risk_score = np.zeros(cox_km_r.shape[0])
for i in range(cox_km_r.shape[0]):
    risk_score[i] = np.sum(np.multiply(coef, gene_ex[i]))
med = np.median(risk_score)
risk_score = np.where(risk_score >= med, 'High risk', 'Low risk')
cox_km_r['Risk'] = risk_score
#cox_km_r.mask(cox_km_r['Risk'] >= med, 'High risk', inplace=True)
cox_km_r.to_csv('fenceng.csv')
gene = cox_km.columns.tolist()

for i in gene_uni:
    if i in l_gene:
        continue
    else:
        cox_km.drop(i, axis=1, inplace=True)
cox_km['Risk'] = risk_score
cox_km.to_csv('fenceng.csv')

mul = pd.read_csv('multivarite_Cox.csv',index_col=0)
coef_ = list(mul.iloc[:, 1])
mul = mul.sort_index()
gene_ex =np.array(cox_km_c.iloc[:, 1:(len(coef) + 1)])
risk_score = np.zeros(cox_km_c.shape[0])
for i in range(cox_km_r.shape[0]):
    risk_score[i] = np.sum(np.multiply(coef, gene_ex[i]))
cox_km_c['Risk'] = risk_score
