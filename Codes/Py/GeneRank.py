# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 19:42:46 2021

@author: 12870
"""
import numpy as np
import pandas as pd
from pandas import Series
from numpy import float32
import math
from numpy import inf
import warnings
import time
import os
from collections import Counter
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier

warnings.filterwarnings("ignore")
expression = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\expr_after_first_selection.csv')
feature = np.array(expression.iloc[:, 1:])
labels = np.array(expression.iloc[:, 0])
expression = np.array(expression)


def get_count_by_counter(l):
    t1 = time.time()
    count = Counter(l)
    t2 = time.time()
    print(t2-t1)
    count_dict = dict(count)
    return count_dict


def get_G():
    tom = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\TOM.csv',index_col=0)
    tom = np.array(tom)
    for i in range(tom.shape[0]):
        tom[i, i] = 0
    add = np.sum(tom)
    number = tom.shape[0]*tom.shape[1]-tom.shape[0]
    mean = add/number
    for i in range(tom.shape[0]):
        for j in range(tom.shape[1]):
            if tom[i, j] >= mean:
                tom[i, j] = 1
            else:
                tom[i, j] = 0
    return tom


def getlog2fc(tag1):
    tumor_1 = train[train[:, 0] == tag1]
    tumor_2 = train[train[:, 0] != tag1]
    tumor_mean1 = np.mean(tumor_1[:, :-1], 0)
    tumor_mean2 = np.mean(tumor_2[:, :-1], 0)
    fc = np.log2(tumor_mean1/tumor_mean2)
    fc[tumor_mean1 == 0] = 0
    fc[tumor_mean2 == 0] = 0
    return fc


def generank(G, fc, d, max_degree=inf):
    fc = np.abs(fc)
    norm_fc = fc / np.max(fc)
    colsum = np.sum(G, 0)
    degrees = np.minimum(max_degree, np.maximum(1, colsum))
    dimG = G.shape[1]
    A = np.eye(dimG)
    D1 = np.zeros((dimG, dimG))
    np.fill_diagonal(D1, 1.0/degrees)
    A = A - d*np.dot(G.T, D1)
    b = (1-d) * norm_fc
    r = np.linalg.solve(A, b)
    return r


def tempranking(r):
    lis = [1]*len(r)
    obj = pd.Series([a / b for a, b in zip(lis, r)])
    order = []
    rr = obj.rank()
    for i in range(len(rr)):
        order.append(rr[i])
    lis2 = [1]*len(order)
    ranks = [a / b for a, b in zip(lis2, order)]
    return ranks


def mkdir(path):
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
    else:
        return False


mkpath = "D:\\project\\TCGA_stomach_analysis\\data"
os.chdir(mkpath)
skf = KFold(n_splits=10, random_state=0, shuffle=True)
d = 0.5
rank = []
G = get_G()
for train_index, test_index in skf.split(feature, labels):
    train = expression[train_index]
    y = labels[train_index]
    num = get_count_by_counter(y)
    tag = train.shape[0]
    # tag = list()
    # for i in range(len(np.unique(y))):
    #   tag.append(num[i+1])
    # min_multi =mcm(tag)
    fc = list()
    r = list()

    for i in range(0, 2):
        fc.append(getlog2fc(i))

    for i in [0, 1]:
        rp = generank(G, fc[i], d, max_degree=inf)
        rp = (rp-min(rp)) / (max(rp)-min(rp))
        coef = num[i] / tag

        r.append(rp*coef)

    rank.append(np.sum(r, 0) / 2)
rank = np.sum(rank, 0) / 10
score = np.array((rank-min(rank))/(max(rank)-min(rank)))
# score=tempranking(rank)
outfile = 'xx.csv'
np.savetxt(outfile, score, delimiter=",")
score = pd.DataFrame(score)
score.to_csv('xx.csv')

max_degree=inf
fc = np.abs(fc[0])
norm_fc = fc / np.max(fc)
colsum = np.sum(G, 0)
degrees = np.minimum(max_degree, np.maximum(1, colsum))
dimG = G.shape[1]
A = np.eye(dimG)
D1 = np.zeros((dimG, dimG))
np.fill_diagonal(D1, 1.0/degrees)
A = A - d*np.dot(G.T, D1)
b = (1-d) * norm_fc.reshape(-1, 1)
r = np.linalg.solve(A, b)


rp = generank(G, fc[i], d, max_degree=inf)
rp = (rp-min(rp)) / (max(rp)-min(rp))
coef = num[i] / tag