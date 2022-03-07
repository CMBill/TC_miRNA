# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 16:50:27 2021

@author: LanHao
"""
from sklearn.linear_model import LogisticRegression
from model_evaluation import plot_roc_curve
from model_evaluation import auc_, mcc_, f1_, precision_, sen_, acc_
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve
from imblearn.over_sampling import SMOTE
import numpy as np
import pandas as pd


# 导入数据 203个肿瘤，24个正常
data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\generank_selected.csv')
X, y = np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0])


clf = LogisticRegression(max_iter=1000)
scores = cross_val_score(clf, X, y, cv=10)
train = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\train.xlsx', header=None)
train = np.array(train)
y_train = train[:, 0]
X_train = train[:, 1:]
test = pd.read_excel('D:\\project\\TCGA_stomach_analysis\\data\\TEST.xlsx', header=None)
test = np.array(test)
y_test = test[:, 0]
X_test = test[:, 1:]
# 网格寻参
tuned_parameters = [{'penalty': ['l2'],
                     'C': [0.001, 0.01, 0.1, 1, 10, 100, 1000],
                     }
                    ]
grid = GridSearchCV(clf,
                    tuned_parameters,
                    cv=10,
                    scoring='neg_log_loss',
                    )
grid.fit(X_train, y_train)
grid.cv_results_
a = grid.best_score_
b = grid.best_params_

# 利用寻找到的的参数建模，评估模型
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    test_size=0.3,
                                                    random_state=0)
model_ = LogisticRegression(penalty='l2', random_state=0, C=1)
model_.fit(X_train, y_train)
prob = model_.predict_proba(X_test)[:, 1]
fper, tper, thresholds = roc_curve(y_test, prob)
plot_roc_curve(fper, tper)
y_pred = model_.predict(X)
auc = auc_(y, y_pred)
mcc = mcc_(y, y_pred)
f1_score = f1_(y, y_pred)
precision = precision_(y, y_pred)
sen = sen_(y, y_pred)
acc = acc_(y, y_pred)
model_eval = pd.DataFrame([dict(auc=auc, mcc=mcc, f1=f1_score,
                                precision=precision,
                                sen=sen, acc=acc)])
model_eval.to_csv('D:\project\TCGA_stomach_analysis\data\download_mirna\milr_eval.csv')


'''
# 得到所有特征的系数
run_times = 100
feature_coef = np.zeros(shape=(100, len(DEGnames)))
for i in range(run_times):
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        random_state=i,
                                                        test_size=0.3)
    clf = LogisticRegression(penalty='l2', random_state=0, C=1)
    model = clf.fit(X_train, y_train)
    importance = model.coef_
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
feature_coef_ = feature_coef_.drop(feature_coef_[(feature_coef_['times'] < 50)].index)
feature_coef_.to_csv(saveNames(2, 0))
'''
