import pandas as pd
import numpy as np
import xgboost as XGB
from sklearn import metrics
from sklearn import svm, preprocessing
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score, confusion_matrix, roc_auc_score, plot_roc_curve, plot_confusion_matrix
from sklearn.metrics import recall_score

path = "E:/wkh/Codes/Projects/TC_miRNA"
# 读取文件
s_count = pd.read_csv(path + '/Data/s_count.csv')
s_rpm = pd.read_csv(path + '/Data/s_rpm.csv')
with open(path + '/Data/degs.txt', 'r') as f:
    degs = []
    for line in f:
        degs.append(line.strip())
with open(path + '/Data/degs.sub.txt', 'r') as f:
    sub = []
    for line in f:
        sub.append(line.strip())
SampleGroup = pd.read_csv(path + '/Data/SampleGroup.csv')
# 数据预处理
SampleGroup.iloc[SampleGroup.iloc[:, 1] == 'cancer', 1] = '1'
SampleGroup.iloc[SampleGroup.iloc[:, 1] == 'normal', 1] = '0'
SampleGroup.index = SampleGroup.iloc[:, 0]
SampleGroup.drop('EntitiesId', axis=1, inplace=True)

rpm_T = s_rpm.T
rpm_T.columns = rpm_T.iloc[0]
rpm_T.drop('NAME', inplace=True)
rpm1 = rpm_T.loc[:, degs]

X_rfe = np.array(rpm1)
y_rfe = np.array(SampleGroup)
# Seek


def svm_seek():
    clf = svm.SVC()
    param_grid = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4], 'C': [1, 10, 100, 1000]},
                  {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]
    grid = GridSearchCV(clf, param_grid, cv=10, scoring='accuracy')
    grid.fit(X_rfe, y_rfe.ravel())
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best


def rf_seek():
    rfc = RFC()
    param_grid = {'n_estimators': [i for i in range(1, 101, 10)],
                  'max_depth': [i for i in range(10, 101, 10)]}
    grid = GridSearchCV(rfc, param_grid, cv=10, scoring='accuracy')
    grid.fit(X_rfe, y_rfe.ravel())
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best


def xgb_seek():
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
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best


svm_p = svm_seek()
rf_p = rf_seek()
xgb_p = xgb_seek()

