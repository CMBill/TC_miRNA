# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 14:26:26 2021

@author: 12870
"""
from model_evaluation import plot_roc_curve
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve,auc
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
import numpy as np
import pandas as pd
import xgboost as xgb

data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\generank_selected.csv')
X, y = np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0])
smo = SMOTE(random_state=42)
X_smo, y_smo = smo.fit_resample(X, y)
X_train, X_test, y_train, y_test = train_test_split(X_smo, y_smo,
                                                    test_size=0.3,
                                                    random_state=0)
clf = xgb.XGBClassifier().fit(X_train, y_train)
model_ = clf
model_.fit(X_train, y_train)
prob = model_.predict_proba(X_test)[:, 1]
fper, tper, thresholds = roc_curve(y_test, prob)
plot_roc_curve(fper, tper)
y_pred = model_.predict(X_smo)
a = auc(fper, tper)
