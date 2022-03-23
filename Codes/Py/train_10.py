# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 10:32:36 2021

@author: 12870
"""
import pandas as pd
import numpy as np
import os
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
from sklearn.metrics import confusion_matrix

os.chdir('D:\\project\\TCGA_stomach_analysis\\data')
vote = pd.read_csv('total_vote.csv',
                   index_col=0)
vote.drop(vote[vote['votes'] < 10].index, inplace=True)
name = vote.index.tolist()
data = pd.read_csv('expr_after_first_selection.csv')
name_ = data.columns.tolist()
del(name_[0])
for i in name_:
    if i in name:
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

rfc = RFC(n_estimators=71, max_depth=10)
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
clf = svm.SVC(kernel='rbf', C=10, gamma=0.001)
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


lr = LogisticRegression()
param_grid = {'penalty': ['l2'],
              'C': [0.001, 0.01, 0.1, 1, 10, 100, 1000]}
grid = GridSearchCV(lr, param_grid, cv=10, scoring='accuracy')
grid.fit(X_rfe, y_rfe)
print(grid.best_score_, grid.best_params_)
lr = LogisticRegression(penalty='l2', C=0.1)
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
    print(cm)
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
knn = KNN(n_neighbors=3, p=5, weights='distance')
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
    print(cm)
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


xgb = XGB.XGBClassifier(n_estimators=6,
                        max_depth=3,
                        min_child_weight=2,
                        gamma=0.1,
                        colsample_bytree=0.8,
                        subsample=0.9,
                        reg_alpha=2,
                        reg_lambda=1)
xgb = XGB.XGBClassifier()
param_grid = {'learning_rate': [0.01, 0.05, 0.07, 0.1, 0.2]}
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
              'max_depth': [i for i in range(10, 101, 10)],}
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
    cm = confusion_matrix(y_true=y_test, y_pred=xgb.predict(X_test))
    print(cm)
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
    cm = confusion_matrix(y_true=y_test, y_pred=dt.predict(X_test))
    print(cm)
    tn, fp, fn, tp = cm.ravel()
    acc.append((tp + tn) / (tp + tn + fp + fn))
    auc.append(roc_auc_score(y_test, dt.predict(X_test)))
    f1.append(f1_score(y_test, dt.predict(X_test)))
    sen.append(tp / (tp + fn))
    spe.append(tn / (tn + fp))
print('acc:{}'.format(np.array(acc).mean()),
      'auc:{}'.format(np.array(auc).mean()),
      'f1:{}'.format(np.array(f1).mean()),
      'spe:{}'.format(np.array(spe).mean()),
      'sen:{}'.format(np.array(sen).mean())
      )
