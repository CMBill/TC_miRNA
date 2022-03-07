# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 14:47:45 2021

@author: 12870
"""

from returnNames import readNames, saveNames
from data_preprocessing import dataProcessing
from readFiles import readFiles
from sklearn.tree import DecisionTreeClassifier as DT
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from model_evaluation import ModelEvaluation
import numpy as np
import pandas as pd

filepath_counts, filepath_names = readNames(eval(input('文件编号：')))
counts_data, DEGnames = readFiles(filepath_counts, filepath_names)
X, y = dataProcessing(counts_data, DEGnames)

clf = DT(criterion='gini',
         random_state=0,
         max_depth=7,
         min_impurity_decrease=0.02631578947,
         min_samples_leaf=1,
         splitter='best')
scores = cross_val_score(clf, X, y, cv=10).mean()

X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state=0,
                                                    test_size=0.3)
param_grid = {'splitter': ['best', 'random'],
              'min_samples_leaf': [i for i in range(1, 50)],
              'min_impurity_decrease': [*np.linspace(0, 0.5, 20)]
              }

grid = GridSearchCV(clf,
                    param_grid,
                    cv=10,
                    )
grid.fit(X_train, y_train)
grid.cv_results_
a = grid.best_score_
b = grid.best_params_


run_times = 100
feature_coef = np.zeros(shape=(100, len(DEGnames)))
for i in range(run_times):
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        random_state=i,
                                                        test_size=0.3)
    clf = DT(criterion='gini',
             random_state=0,
             max_depth=7,
             min_impurity_decrease=0.02631578947,
             min_samples_leaf=1,
             splitter='best')
    model = clf.fit(X_train, y_train)
    importance = model.feature_importances_
    feature_coef[i] = importance

feature_coef = feature_coef.T
feature_coef = np.concatenate((feature_coef,
                               np.mean(feature_coef, axis=1).reshape(-1, 1)),
                              axis=1)
runtime = [i for i in range(1, 101)] + ['mean']
feature_coef_ = pd.DataFrame(feature_coef,
                             index=DEGnames,
                             columns=runtime)

# 对得到的特征系数矩阵进行筛选
threshold = np.mean(feature_coef)
counts = 0
importances = np.sum(feature_coef[:, :-1] > threshold, axis=1)
importances_ = pd.DataFrame(importances, index=DEGnames, columns=['times'])
feature_coef_ = pd.concat([feature_coef_, importances_], axis=1)
feature_coef_ = feature_coef_.drop(feature_coef_[(feature_coef_['times'] <= 0)].index)
feature_coef_.to_csv(saveNames(2, 5))
