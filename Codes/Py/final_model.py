# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 18:45:42 2021

@author: 12870
"""


import pandas as pd
import numpy as np
import os
import re
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import RFE
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier as RFC
from collections import Counter
from sklearn import svm
import xgboost as XGB
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier as DT
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier as DT
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score, recall_score, f1_score
from imblearn.over_sampling import SMOTE
from sklearn.metrics import confusion_matrix

data = pd.read_csv('after_rfe.csv')
X_rfe = np.array(data.iloc[:, 1:])
y_rfe = np.array(data.iloc[:, 0])
xgb = XGB.XGBClassifier(n_estimators=24,
                        max_depth=4,
                        min_child_weight=1,
                        colsample_bytree=0.7,
                        subsample=0.7,
                        reg_alpha=0.05,
                        reg_lambda=1)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    xgb.fit(X_train, y_train)
    cm = confusion_matrix(y_true=y_test, y_pred=xgb.predict(X_test))
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, xgb.predict(X_test)))
    f1.append(f1_score(y_test, xgb.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )


lr = LogisticRegression(penalty='l2', C=10, max_iter=1000)
lr = LogisticRegression()
param_grid = {'penalty': ['l2'],
              'C': [0.001, 0.01, 0.1, 1, 10, 100, 1000]}
grid = GridSearchCV(lr, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    lr.fit(X_train, y_train)
    cm = confusion_matrix(y_true=y_test, y_pred=lr.predict(X_test))
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, lr.predict(X_test)))
    f1.append(f1_score(y_test, lr.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

data_ = pd.read_csv('de_val.csv', index_col=0)
X, y = np.array(data_.iloc[:, 1:]), np.array(data_.iloc[:, 0])
skf = StratifiedKFold(n_splits=5)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X, y):
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    lr.fit(X_train, y_train)
    cm = confusion_matrix(y_true=y_test, y_pred=lr.predict(X_test))
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, lr.predict(X_test)))
    f1.append(f1_score(y_test, lr.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )
lr.fit(X, y)
cm = confusion_matrix(y_true=y, y_pred=lr.predict(X))
tn, fp, fn, tp = cm.ravel()
acc.append((tp + tn) / (tp + tn + fp + fn))
auc.append(roc_auc_score(y, lr.predict(X)))
f1.append(f1_score(y, lr.predict(X)))
sen.append(tp / (tp + fn))
spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )


knn = KNN(n_neighbors=3, weights='uniform')
param_grid = [
    {
        'weights': ['uniform'],
        'n_neighbors': [i for i in range(1, 11)]
    },
    {
        'weights': ['distance'],
        'n_neighbors': [i for i in range(1, 11)],
        'p': [i for i in range(1, 6)]
    }
]
grid = GridSearchCV(knn, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)

skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    knn.fit(X_train, y_train)
    cm = confusion_matrix(y_true=y_test, y_pred=knn.predict(X_test))
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, knn.predict(X_test)))
    f1.append(f1_score(y_test, knn.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

skf = StratifiedKFold(n_splits=5)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X, y):
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    knn.fit(X_train, y_train)
    cm = confusion_matrix(y_true=y_test, y_pred=knn.predict(X_test))
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, knn.predict(X_test)))
    f1.append(f1_score(y_test, knn.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )
miRNA = pd.read_csv('miRNA.csv', index_col=0)
val_name = miRNA.columns.tolist()
sample_l = []
for i in val_name:
    if re.match(r'^TCGA\.\w\w\.\w\w\w\w\.01A', i):
        i = 1
        sample_l.append(i)
    else:
        i = 0
        sample_l.append(i)
f = pd.read_csv('D:\project\cancer subtype\data\cnv\TCGA.STAD.sampleMap_SNP6_nocnv_genomicSegment\SNP6_nocnv_genomicSegment.tsv')
