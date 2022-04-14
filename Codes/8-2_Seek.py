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
    grid = GridSearchCV(clf, param_grid, cv=5, scoring='accuracy')
    grid.fit(X_rfe, y_rfe.ravel())
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9544927536231885, {'C': 1, 'kernel': 'linear'}]


def rf_seek():
    rfc = RandomForestClassifier()
    param_grid = {'n_estimators': [i for i in range(1, 101, 10)],
                  'max_depth': [i for i in range(10, 101, 10)]}
    grid = GridSearchCV(rfc, param_grid, cv=5, scoring='accuracy')
    grid.fit(X_rfe, y_rfe.ravel())
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9895194508009155, {'max_depth': 70, 'n_estimators': 61}]


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
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9842868039664378, {'max_depth': 3, 'min_child_weight': 1}]


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
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9842868039664378, {'gamma': 0.0}]


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
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9877803203661328, {'colsample_bytree': 0.6, 'subsample': 0.6}]


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
                                 subsample=0.6,
                                 colsample_bytree=0.6)
    param_test1 = {'reg_alpha': [1e-5, 1e-2, 0.1, 1, 100]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9877955758962624, {'reg_alpha': 0.01}]


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
                                 subsample=0.6,
                                 colsample_bytree=0.6,
                                 reg_alpha=0.01)
    param_test1 = {'reg_lambda': [1e-5, 1e-2, 0.1, 1, 100]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9877955758962624, {'reg_lambda': 1}]


def xgb_seek6():
    xgb1 = xgboost.XGBClassifier(n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=3,
                                 min_child_weight=1,
                                 gamma=0,
                                 subsample=0.6,
                                 colsample_bytree=0.6,
                                 reg_alpha=0.01,
                                 reg_lambda=1)
    param_test1 = {'learning_rate': [0.01, 0.015, 0.025, 0.05, 0.1]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9877955758962624, {'learning_rate': 0.1}]


def xgb_seek7():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=3,
                                 min_child_weight=1,
                                 gamma=0,
                                 subsample=0.6,
                                 colsample_bytree=0.6,
                                 reg_alpha=0.01,
                                 reg_lambda=1)
    param_test1 = {'n_estimators': [i for i in range(10, 500, 10)]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(X_rfe, y_rfe)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.9895347063310449, {'n_estimators': 170}]


svm_p = svm_seek()
rf_p = rf_seek()
xgb_p1 = xgb_seek1()
xgb_p2 = xgb_seek2()
xgb_p3 = xgb_seek3()
xgb_p4 = xgb_seek4()
xgb_p5 = xgb_seek5()
xgb_p6 = xgb_seek6()
xgb_p7 = xgb_seek7()

print('svm: ' + str(svm_p) + '\n'
      + 'rf: ' + str(rf_p) + '\n'
      + 'xgb1: ' + str(xgb_p1) + '\n'
      + 'xgb2: ' + str(xgb_p2) + '\n'
      + 'xgb3: ' + str(xgb_p3) + '\n'
      + 'xgb4: ' + str(xgb_p4) + '\n'
      + 'xgb5: ' + str(xgb_p5) + '\n'
      + 'xgb6: ' + str(xgb_p6) + '\n'
      + 'xgb7: ' + str(xgb_p7) + '\n')
'''
svm: [0.9544927536231885, {'C': 1, 'kernel': 'linear'}]
rf: [0.9895194508009155, {'max_depth': 70, 'n_estimators': 61}]
xgb1: [0.9842868039664378, {'max_depth': 3, 'min_child_weight': 1}]
xgb2: [0.9842868039664378, {'gamma': 0.0}]
xgb3: [0.9877803203661328, {'colsample_bytree': 0.6, 'subsample': 0.6}]
xgb4: [0.9877955758962624, {'reg_alpha': 0.01}]
xgb5: [0.9877955758962624, {'reg_lambda': 1}]
xgb6: [0.9877955758962624, {'learning_rate': 0.1}]
xgb7: [0.9895347063310449, {'n_estimators': 170}]
'''