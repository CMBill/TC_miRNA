import pandas as pd
import numpy as np
import xgboost
from sklearn import svm
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
import xgboost as XGB
from sklearn.ensemble import RandomForestClassifier as RFC
import Codes.fit_fc as fc

path = "."
# 读取文件
s_count = pd.read_csv(path + '/Data/s_count.csv')
s_rpm = pd.read_csv(path + '/Data/s_rpm.csv')
with open(path + '/Data/degs.sub.txt', 'r') as f:
    sub = []
    for line in f:
        sub.append(line.strip())
sub_group = pd.read_csv(path + '/Data/sub_group.csv')
# 数据预处理
s_group_ptc = sub_group.loc[sub_group.loc[:, 'Group'] == 'PTC', ]
s_group_fvptc = sub_group.loc[sub_group.loc[:, 'Group'] == 'FV-PTC', ]
s_group_n = sub_group.loc[sub_group.loc[:, 'Group'] == 'normal', ]
s_group = pd.concat([s_group_ptc, s_group_fvptc])

s_group.iloc[s_group.iloc[:, 1] == 'PTC', 1] = '1'
s_group.iloc[s_group.iloc[:, 1] == 'FV-PTC', 1] = '0'
s_group.index = s_group.iloc[:, 0]
s_group.drop('EntitiesId', axis=1, inplace=True)

rpm_T = s_rpm.T
rpm_T.columns = rpm_T.iloc[0]
rpm_T.drop('NAME', inplace=True)
rpm1 = rpm_T.loc[:, sub]

rpm2 = rpm1.loc[s_group.index]

X_rfe = np.array(rpm2)
y_rfe = np.array(s_group)

x_train, x_test, y_train, y_test = train_test_split(rpm1, s_group, random_state=100, stratify=s_group,
                                                     train_size=0.6, test_size=0.4)

def svm_seek():
    clf = svm.SVC()
    # param_grid = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4], 'C': [1, 10, 100, 1000]}]
    param_grid = [{'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]
    grid = GridSearchCV(clf, param_grid, cv=5, scoring='accuracy')
    grid.fit(x_train, y_train.ravel())
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best


def rf_seek():
    rfc = RandomForestClassifier()
    param_grid = {'n_estimators': [i for i in range(1, 101, 10)],
                  'max_depth': [i for i in range(10, 101, 10)]}
    grid = GridSearchCV(rfc, param_grid, cv=5, scoring='accuracy')
    grid.fit(x_train, y_train.ravel())
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8525966598032486, {'max_depth': 10, 'n_estimators': 71}]

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
    grid.fit(x_train, y_train)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8398078242964997, {'max_depth': 5, 'min_child_weight': 3}]


def xgb_seek2():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 subsample=0.8,
                                 colsample_bytree=0.8,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=5,
                                 min_child_weight=3)
    param_test1 = {'gamma': [i / 10.0 for i in range(0, 5)]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(x_train, y_train)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8398078242964997, {'gamma': 0.0}]


def xgb_seek3():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=5,
                                 min_child_weight=3,
                                 gamma=0.0)
    param_test1 = {'subsample': [i / 10.0 for i in range(6, 10)],
                   'colsample_bytree': [i / 10.0 for i in range(6, 10)]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(x_train, y_train)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8461450469000228, {'colsample_bytree': 0.6, 'subsample': 0.6}]


def xgb_seek4():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=5,
                                 min_child_weight=3,
                                 gamma=0,
                                 subsample=0.6,
                                 colsample_bytree=0.6)
    param_test1 = {'reg_alpha': [1e-5, 1e-2, 0.1, 1, 100]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(x_train, y_train)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8461679249599634, {'reg_alpha': 0.1}]


def xgb_seek5():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=5,
                                 min_child_weight=3,
                                 gamma=0,
                                 subsample=0.6,
                                 colsample_bytree=0.6,
                                 reg_alpha=0.1)
    param_test1 = {'reg_lambda': [1e-5, 1e-2, 0.1, 1, 100]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(x_train, y_train)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8461679249599634, {'reg_lambda': 1}]


def xgb_seek6():
    xgb1 = xgboost.XGBClassifier(n_estimators=140,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=5,
                                 min_child_weight=3,
                                 gamma=0,
                                 subsample=0.6,
                                 colsample_bytree=0.6,
                                 reg_alpha=0.1,
                                 reg_lambda=1)
    param_test1 = {'learning_rate': [0.01, 0.015, 0.025, 0.05, 0.1]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(x_train, y_train)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8461679249599634, {'learning_rate': 0.1}]


def xgb_seek7():
    xgb1 = xgboost.XGBClassifier(learning_rate=0.1,
                                 objective='binary:logistic',
                                 nthread=4,
                                 scale_pos_weight=1,
                                 seed=27,
                                 max_depth=5,
                                 min_child_weight=3,
                                 gamma=0,
                                 subsample=0.6,
                                 colsample_bytree=0.6,
                                 reg_alpha=0.1,
                                 reg_lambda=1)
    param_test1 = {'n_estimators': [i for i in range(10, 500, 10)]}
    grid = GridSearchCV(xgb1, param_test1, cv=5, scoring='accuracy', n_jobs=4)
    grid.fit(x_train, y_train)
    param_best = [grid.best_score_, grid.best_params_]
    print(param_best)
    return param_best
# [0.8547014413177763, {'n_estimators': 20}]


svm_p = svm_seek()
# [0.7756706753006475, {'C': 1, 'gamma': 0.001, 'kernel': 'rbf'}]
rf_p = rf_seek()
# [0.8503237742830713, {'max_depth': 80, 'n_estimators': 61}]
# [0.8567067530064755, {'max_depth': 20, 'n_estimators': 61}]
# [0.8524514338575393, {'max_depth': 100, 'n_estimators': 51}]
# [0.8523589269195189, {'max_depth': 60, 'n_estimators': 61}]
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
      + 'xgb5: ' + str(xgb_p5) + '\n')

#fit
print('svm:')
clf = svm.SVC(kernel='linear', C=1, probability=True)  # kernel='linear', C=1
fc.fit_c(clf, x_test, y_test, 5)
print('rf:')
rfc = RFC(n_estimators=61, max_depth=20)
fc.fit_c(rfc, x_test, y_test, 5)
print('xgb:')
xgb = XGB.XGBClassifier(objective='binary:logistic',
                        nthread=4,
                        learning_rate=0.1,
                        n_estimators=20,
                        max_depth=5,
                        min_child_weight=3,
                        gamma=0,
                        colsample_bytree=0.6,
                        subsample=0.6,
                        reg_alpha=0.1,
                        reg_lambda=1)
fc.fit_c(xgb, x_test, y_test, 5)
