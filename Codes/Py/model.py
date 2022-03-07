# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:00:06 2021

@author: 12870
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier as DT
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from collections import Counter
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score, recall_score, f1_score

os.chdir('D:\\project\\TCGA_stomach_analysis\\data')
data = pd.read_excel('generank_selected.xlsx')
X, y = np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0])
train = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\train.xlsx', header=None)
train = np.array(train)
y_train = train[:, 0]
X_train = train[:, 1:]
test = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\TEST.xlsx', header=None)
test = np.array(test)
y_test = test[:, 0]
X_test = test[:, 1:]
# 对每一个算法先调参
SVM = svm.SVC(kernel='rbf', gamma=0.0001, C=1000).fit(X_train, y_train)
score = cross_val_score(SVM, X, y, cv=10, scoring='roc_auc').mean()
cross_val_score(SVM, X, y, cv=10).mean()
tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
                     'C': [1, 10, 100, 1000]},
                    {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]
grid = GridSearchCV(SVM,
                    tuned_parameters,
                    cv=10,
                    scoring='roc_auc',
                    )
grid.fit(X_train, y_train)
grid.cv_results_
a = grid.best_score_
b = grid.best_params_
print('测试集acc：', accuracy_score(y_test, SVM.predict(X_test)))
print('测试集SPE：',precision_score(y_test, SVM.predict(X_test)))
print('测试集AUC：',roc_auc_score(y_test, SVM.predict(X_test)))
print('测试集f1:',f1_score(y_test, SVM.predict(X_test)))
print('测试集sen:',recall_score(y_test, SVM.predict(X_test)))

rfc = RFC(n_estimators=75, max_depth=12, min_samples_split=2).fit(X_train,
                                                                     y_train)
print('测试集acc：', accuracy_score(y_test, rfc.predict(X_test)))
print('测试集SPE：',precision_score(y_test, rfc.predict(X_test)))
print('测试集AUC：',roc_auc_score(y_test, rfc.predict(X_test)))
print('测试集f1:',f1_score(y_test, rfc.predict(X_test)))
print('测试集sen:',recall_score(y_test, rfc.predict(X_test)))

lr = LogisticRegression(penalty='l2', random_state=0, C=1).fit(X_train, y_train)
score = cross_val_score(lr, X, y, cv=10, scoring='roc_auc').mean()
print('测试集acc：', accuracy_score(y_test, lr.predict(X_test)))
print('测试集SPE：',precision_score(y_test, lr.predict(X_test)))
print('测试集AUC：',roc_auc_score(y_test, lr.predict(X_test)))
print('测试集f1:',f1_score(y_test, lr.predict(X_test)))
print('测试集sen:',recall_score(y_test, lr.predict(X_test)))

clf = xgb.XGBClassifier().fit(X_train, y_train)
print('测试集acc：', accuracy_score(y_test, clf.predict(X_test)))
print('测试集SPE：',precision_score(y_test, clf.predict(X_test)))
print('测试集AUC：',roc_auc_score(y_test, clf.predict(X_test)))
print('测试集f1:',f1_score(y_test, clf.predict(X_test)))
print('测试集sen:',recall_score(y_test, clf.predict(X_test)))

dt = DT(criterion='gini',min_samples_leaf=2)
dt = DT(criterion='gini',
         random_state=0,
         max_depth=30,
         min_impurity_decrease=0.2,
         min_samples_leaf=2,
         splitter='best').fit(X_train, y_train)
param = [{'criterion':['gini'],'max_depth':[i for i in range(31)],'min_samples_leaf':[2,3,5,10],'min_impurity_decrease':[0.1,0.2,0.5]},
        ]
grid = GridSearchCV(dt, param_grid=param,cv=10)
grid.fit(X_train,y_train)
print('最优分类器:',grid.best_params_,'最优分数:', grid.best_score_)

print('测试集acc：', accuracy_score(y_test, dt.predict(X_test)))
print('测试集SPE：',precision_score(y_test, dt.predict(X_test)))
print('测试集AUC：',roc_auc_score(y_test, dt.predict(X_test)))
print('测试集f1:',f1_score(y_test, dt.predict(X_test)))
print('测试集sen:',recall_score(y_test, dt.predict(X_test)))

knn = KNeighborsClassifier(n_neighbors=7, p=3, weights='distance').fit(X_train,
                                                                       y_train)
print('测试集acc：', accuracy_score(y_test, knn.predict(X_test)))
print('测试集SPE：',precision_score(y_test, knn.predict(X_test)))
print('测试集AUC：',roc_auc_score(y_test, knn.predict(X_test)))
print('测试集f1:',f1_score(y_test, knn.predict(X_test)))
print('测试集sen:',recall_score(y_test, knn.predict(X_test)))

data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\final.csv')
X, y = np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0])
train = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\finalTrain.xlsx',header=0)
train = np.array(train)
X_train = train[:, 1:]
y_train = train[:, 0]
test = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\finalTest.xlsx',header=0)
test = np.array(test)
X_test = test[:, 1:]
y_test = test[:, 0]
lr = LogisticRegression().fit(X_train, y_train)
print('测试集acc：', accuracy_score(y_test, lr.predict(X_test)))
print('测试集SPE：',precision_score(y_test, lr.predict(X_test)))
print('测试集AUC：',roc_auc_score(y_test, lr.predict(X_test)))
print('测试集f1:',f1_score(y_test, lr.predict(X_test)))
print('测试集sen:',recall_score(y_test, lr.predict(X_test)))
