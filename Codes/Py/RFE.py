# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 11:26:25 2021

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
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score, recall_score, f1_score
from imblearn.over_sampling import SMOTE
from sklearn.metrics import confusion_matrix

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

X, y = np.array(rna.iloc[:, 1:]), np.array(rna.iloc[:, 0])


iterations = 10
gene_name = rna.columns.tolist()[1:]
rfc = RFC()
clf = svm.SVC(kernel='linear')
gene_rank_rf = pd.DataFrame(index=gene_name)
gene_rank_svm = pd.DataFrame(index=gene_name)
for i in range(1, iterations + 1):
    rfecv = RFECV(estimator=rfc,
                  cv=StratifiedKFold(n_splits=10, shuffle=True),
                  scoring='accuracy')
    rfecv.fit(X, y)
    rank = rfecv.ranking_
    gene_rank_rf['round_{}'.format(i)] = rank
    rfecv_ = RFECV(estimator=clf,
                   cv=StratifiedKFold(n_splits=10, shuffle=True),
                   scoring='accuracy')
    rfecv_.fit(X, y)
    rank_ = rfecv_.ranking_
    gene_rank_svm['round_{}'.format(i)] = rank_

gene_rank_rf['votes'] = (gene_rank_rf == 1).sum(axis=1)
gene_rank_rf = gene_rank_rf.T
gene_rank_rf['ss'] = (gene_rank_rf == 1).sum(axis=1)
#gene_rank_rf = gene_rank_rf[~gene_rank_rf.isin([0])].dropna(axis=1)
a = dict(gene_rank_rf.loc['votes'].value_counts())

gene_rank_svm['votes'] = (gene_rank_svm == 1).sum(axis=1)
gene_rank_svm = gene_rank_svm.T
gene_rank_svm['ss'] = (gene_rank_svm == 1).sum(axis=1)
#gene_rank_svm = gene_rank_svm[~gene_rank_svm.isin([0])].dropna(axis=1)
b = dict(gene_rank_svm.loc['votes'].value_counts())
total_vote = pd.DataFrame(gene_rank_rf.loc['votes'] + gene_rank_svm.loc['votes'])
c = dict(total_vote.iloc[:, 0].value_counts())
vote_count = pd.DataFrame(columns=['counts', 'votes'])
vote_count['counts'] = list(c.values())
vote_count['votes'] = list(c.keys())
rfe_name = total_vote.index.tolist()[:-1]
for i in gene_name:
    if i in rfe_name:
        continue
    else:
        rna.drop(i, axis=1, inplace=True)

data = pd.read_csv('rfe_selection.csv', index_col=0)
data = rna.copy()
X_rfe = np.array(data.iloc[:, 1:])
y_rfe = np.array(data.iloc[:, 0])
smo = SMOTE(random_state=0)
X_smo, y_smo = smo.fit_resample(X_rfe, y_rfe)

total_vote.drop(['ss'], axis=0, inplace=True)
vote = total_vote.copy()
data.to_csv('after_rfe.csv')


vote.drop(vote[vote.votes < 11].index, inplace=True)
vote.to_csv('enrich.csv')
sel_name = vote.index.tolist()
for i in data.columns.tolist()[1:]:
    if i in sel_name:
        continue
    else:
        data.drop(i, axis=1, inplace=True)
X_rfe = np.array(data.iloc[:, 1:])
y_rfe = np.array(data.iloc[:, 0])
rfc = RFC()
param_grid = {'n_estimators': [i for i in range(1, 101, 10)],
              'max_depth': [i for i in range(10, 101, 10)],
              }
grid = GridSearchCV(rfc, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_) 

rfc = RFC(n_estimators=1, max_depth=30)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    rfc.fit(X_train, y_train)
    cm = confusion_matrix(y_true=y_test, y_pred=rfc.predict(X_test))
    print(cm)
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, rfc.predict(X_test)))
    f1.append(f1_score(y_test, rfc.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

clf = svm.SVC()
param_grid = [{'kernel': ['rbf'],
               'gamma': [1e-3, 1e-4],
               'C': [1, 10, 100, 1000]
               },
              {'kernel': ['linear'],
               'C': [1, 10, 100, 1000]
               }
              ]
grid = GridSearchCV(clf, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)
clf = svm.SVC(kernel='rbf', C=100, gamma=0.001)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    clf.fit(X_train, y_train)
    cm = confusion_matrix(y_true=y_test, y_pred=clf.predict(X_test))
    print(cm)
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, clf.predict(X_test)))
    f1.append(f1_score(y_test, clf.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

lr = LogisticRegression(penalty='l2', C=0.1)
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
    acc.append(accuracy_score(y_test, lr.predict(X_test)))
    auc.append(roc_auc_score(y_test, lr.predict(X_test)))
    f1.append(f1_score(y_test, lr.predict(X_test)))
    sen.append(recall_score(y_test, lr.predict(X_test)))
    spe.append(precision_score(y_test, lr.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

knn = KNN()
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
knn = KNN(n_neighbors=2, weights='uniform')
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
    acc.append(accuracy_score(y_test, knn.predict(X_test)))
    auc.append(roc_auc_score(y_test, knn.predict(X_test)))
    f1.append(f1_score(y_test, knn.predict(X_test)))
    sen.append(recall_score(y_test, knn.predict(X_test)))
    spe.append(precision_score(y_test, knn.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )


xgb = XGB.XGBClassifier(n_estimators=6,
                        max_depth=3,
                        min_child_weight=2,
                        colsample_bytree=0.6,
                        subsample=0.8,
                        reg_alpha=0.05,
                        reg_lambda=0.1,
                        learning_rate=0.2)
param_grid = {'learning_rate': [0.01, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5]}
grid = GridSearchCV(xgb, param_grid, cv=10, scoring='accuracy')
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
    xgb.fit(X_train, y_train)
    acc.append(accuracy_score(y_test, xgb.predict(X_test)))
    auc.append(roc_auc_score(y_test, xgb.predict(X_test)))
    f1.append(f1_score(y_test, xgb.predict(X_test)))
    sen.append(recall_score(y_test, xgb.predict(X_test)))
    spe.append(precision_score(y_test, xgb.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )


dt = DT()
param_grid = {'criterion': ['gini', 'entropy'],
              'max_depth': [i for i in range(1, 31, 5)],
              'min_samples_leaf': [2, 3, 5, 10],
              'min_impurity_decrease': [0.1, 0.2, 0.5]}
grid = GridSearchCV(dt, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)
dt = DT(criterion='entropy',
        max_depth=1,
        min_impurity_decrease=0.1,
        min_samples_leaf=2)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    dt.fit(X_train, y_train)
    acc.append(accuracy_score(y_test, dt.predict(X_test)))
    auc.append(roc_auc_score(y_test, dt.predict(X_test)))
    f1.append(f1_score(y_test, dt.predict(X_test)))
    sen.append(recall_score(y_test, dt.predict(X_test)))
    spe.append(precision_score(y_test, dt.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )


data = pd.read_csv('rfe_selection.csv', index_col=0)
data.drop('hsa.mir.196a.1', axis=1, inplace=True)
X_rfe = np.array(data.iloc[:, 1:])
y_rfe = np.array(data.iloc[:, 0])

rfc = RFC()
param_grid = {'n_estimators': [i for i in range(1, 101, 10)],
              'max_depth': [i for i in range(10, 101, 10)],
              }
grid = GridSearchCV(rfc, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)

rfc = RFC(n_estimators=51, max_depth=10)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    rfc.fit(X_train, y_train)
    acc.append(accuracy_score(y_test, rfc.predict(X_test)))
    auc.append(roc_auc_score(y_test, rfc.predict(X_test)))
    f1.append(f1_score(y_test, rfc.predict(X_test)))
    sen.append(recall_score(y_test, rfc.predict(X_test)))
    spe.append(precision_score(y_test, rfc.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

clf = svm.SVC()
param_grid = [{'kernel': ['rbf'],
               'gamma': [1e-3, 1e-4],
               'C': [1, 10, 100, 1000]
               },
              {'kernel': ['linear'],
               'C': [1, 10, 100, 1000]
               }
              ]
grid = GridSearchCV(clf, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)
clf = svm.SVC(kernel='rbf', C=100, gamma=0.001)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    clf.fit(X_train, y_train)
    acc.append(accuracy_score(y_test, clf.predict(X_test)))
    auc.append(roc_auc_score(y_test, clf.predict(X_test)))
    f1.append(f1_score(y_test, clf.predict(X_test)))
    sen.append(recall_score(y_test, clf.predict(X_test)))
    spe.append(precision_score(y_test, clf.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

lr = LogisticRegression(penalty='l2', C=0.1)
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
    acc.append(accuracy_score(y_test, lr.predict(X_test)))
    auc.append(roc_auc_score(y_test, lr.predict(X_test)))
    f1.append(f1_score(y_test, lr.predict(X_test)))
    sen.append(recall_score(y_test, lr.predict(X_test)))
    spe.append(precision_score(y_test, lr.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )

knn = KNN()
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
knn = KNN(n_neighbors=1, p=4, weights='distance')
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
    acc.append(accuracy_score(y_test, knn.predict(X_test)))
    auc.append(roc_auc_score(y_test, knn.predict(X_test)))
    f1.append(f1_score(y_test, knn.predict(X_test)))
    sen.append(recall_score(y_test, knn.predict(X_test)))
    spe.append(precision_score(y_test, knn.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )


xgb = XGB.XGBClassifier(n_estimators=6,
                        max_depth=3,
                        min_child_weight=2,
                        colsample_bytree=0.6,
                        subsample=0.8,
                        reg_alpha=0.05,
                        reg_lambda=0.1,
                        learning_rate=0.2)
param_grid = {'learning_rate': [0.01, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5]}
grid = GridSearchCV(xgb, param_grid, cv=10, scoring='accuracy')
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
    xgb.fit(X_train, y_train)
    acc.append(accuracy_score(y_test, xgb.predict(X_test)))
    auc.append(roc_auc_score(y_test, xgb.predict(X_test)))
    f1.append(f1_score(y_test, xgb.predict(X_test)))
    sen.append(recall_score(y_test, xgb.predict(X_test)))
    spe.append(precision_score(y_test, xgb.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )


dt = DT()
param_grid = {'criterion': ['gini', 'entropy'],
              'max_depth': [i for i in range(1, 31, 5)],
              'min_samples_leaf': [2, 3, 5, 10],
              'min_impurity_decrease': [0.1, 0.2, 0.5]}
grid = GridSearchCV(dt, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)
dt = DT(criterion='entropy',
        max_depth=1,
        min_impurity_decrease=0.1,
        min_samples_leaf=2)
skf = StratifiedKFold(n_splits=10)
acc = []
auc = []
f1 = []
spe = []
sen = []
for train_index, test_index in skf.split(X_rfe, y_rfe):
    X_train, X_test = X_rfe[train_index], X_rfe[test_index]
    y_train, y_test = y_rfe[train_index], y_rfe[test_index]
    dt.fit(X_train, y_train)
    acc.append(accuracy_score(y_test, dt.predict(X_test)))
    auc.append(roc_auc_score(y_test, dt.predict(X_test)))
    f1.append(f1_score(y_test, dt.predict(X_test)))
    sen.append(recall_score(y_test, dt.predict(X_test)))
    spe.append(precision_score(y_test, dt.predict(X_test)))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )
