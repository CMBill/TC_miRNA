# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 09:38:28 2021

@author: 12870
"""
from returnNames import readNames, saveNames
from data_preprocessing import dataProcessing
from readFiles import readFiles
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from model_evaluation import ModelEvaluation
import numpy as np
import pandas as pd

filepath_counts, filepath_names = readNames(eval(input('文件编号：')))
counts_data, DEGnames = readFiles(filepath_counts, filepath_names)
X, y = dataProcessing(counts_data, DEGnames)

clf = KNeighborsClassifier(n_neighbors=5).fit(X, y)
scores = cross_val_score(clf, X, y, cv=10)

X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state=0,
                                                    test_size=0.3)
param_grid = [
    {'weights': ['uniform'],
     'n_neighbors':[i for i in range(1, 11)]
     },
    {'weights': ['distance'],
     'n_neighbors':[i for i in range(1, 11)],
     'p':[i for i in range(1, 11)]
     }
    ]
knn_clf = KNeighborsClassifier()
grid = GridSearchCV(knn_clf, param_grid)
grid.fit(X_train, y_train)
a = grid.best_estimator_
b = grid.best_params_
c = grid.best_index_
d = grid.best_score_
