import pandas as pd
import numpy as np
import xgboost
from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

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
    rfc = RandomForestClassifier()
    param_grid = {'n_estimators': [i for i in range(1, 101, 10)],
                  'max_depth': [i for i in range(10, 101, 10)]}
    grid = GridSearchCV(rfc, param_grid, cv=10, scoring='accuracy')
    grid.fit(X_rfe, y_rfe.ravel())
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best


def xgb_seek1():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 gamma=0,
                                 subsample=0.8,
                                 colsample_bytree=0.8,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27)
    param_test1 = {'max_depth': [i for i in range(3, 10, 2)],
                   'min_child_weight': [i for i in range(1, 6, 2)]}
    grid = GridSearchCV(xgb1, param_test1, cv=10, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9895039322444041, {'max_depth': 3, 'min_child_weight': 1}]


def xgb_seek2():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 subsample=0.8,
                                 colsample_bytree=0.8,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=3,
                                 min_child_weight=1)
    param_test1 = {'gamma': [i / 10.0 for i in range(0, 5)]}
    grid = GridSearchCV(xgb1, param_test1, cv=10, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9895039322444041, {'gamma': 0.0}]


def xgb_seek3():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=3,
                                 min_child_weight=1,
                                 gamma=0)
    param_test1 = {'subsample': [i / 10.0 for i in range(6, 10)],
                   'colsample_bytree': [i / 10.0 for i in range(6, 10)]}
    grid = GridSearchCV(xgb1, param_test1, cv=10, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9895039322444041, {'colsample_bytree': 0.8, 'subsample': 0.8}]


def xgb_seek4():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=3,
                                 min_child_weight=1,
                                 gamma=0,
                                 subsample=0.8,
                                 colsample_bytree=0.8)
    param_test1 = {'reg_alpha': [1e-5, 1e-2, 0.1, 1, 100]}
    grid = GridSearchCV(xgb1, param_test1, cv=10, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9895039322444041, {'reg_alpha': 1e-05}]


def xgb_seek5():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=3,
                                 min_child_weight=1,
                                 gamma=0,
                                 subsample=0.8,
                                 colsample_bytree=0.8,
                                 reg_alpha=1e-5)
    param_test1 = {'reg_lambda': [1e-5, 1e-2, 0.1, 1, 100]}
    grid = GridSearchCV(xgb1, param_test1, cv=10, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9895039322444041, {'reg_lambda': 1}]


svm_p = svm_seek()
rf_p = rf_seek()
xgb_p1 = xgb_seek1()
xgb_p2 = xgb_seek2()
xgb_p3 = xgb_seek3()
xgb_p4 = xgb_seek4()
xgb_p5 = xgb_seek5()
