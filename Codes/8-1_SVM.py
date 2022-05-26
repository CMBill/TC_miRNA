import pandas as pd
import numpy as np
import xgboost as XGB
from sklearn import metrics
from sklearn import svm
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score, confusion_matrix, roc_auc_score, plot_roc_curve, plot_confusion_matrix
from sklearn.metrics import recall_score

# 读取文件
path = "E:/wkh/Codes/Projects/TC_miRNA"
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
SampleGroup.iloc[SampleGroup.iloc[:, 1] == 'cancer', 1] = '1'
SampleGroup.iloc[SampleGroup.iloc[:, 1] == 'normal', 1] = '0'

SampleGroup.index = SampleGroup.iloc[:, 0]
SampleGroup.drop('EntitiesId', axis=1, inplace=True)
# 选出差异表达的miRNA
rpm_T = s_rpm.T
rpm_T.columns = rpm_T.iloc[0]
rpm_T.drop('NAME', inplace=True)
rpm1 = rpm_T.loc[:, degs]

X_rfe = np.array(rpm1)
y_rfe = np.array(SampleGroup)
# 划分训练测试集
# x_train, x_test, y_train, y_test = train_test_split(rpm1, SampleGroup, random_state=100, stratify=SampleGroup,
#                                                     train_size=0.6, test_size=0.4)
# SVM
def svm_cv():
    clf = svm.SVC()
    param_grid = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4], 'C': [1, 10, 100, 1000]},
                  {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]
    grid = GridSearchCV(clf, param_grid, cv=10, scoring='accuracy')
    grid.fit(rpm1, SampleGroup.values.ravel())
    print(grid.best_score_, grid.best_params_)
    clf = svm.SVC(kernel='linear', C=1)
    skf = StratifiedKFold(n_splits=10)
    acc = []
    auc = []
    spe = []
    sen = []
    for train_index, test_index in skf.split(rpm1, SampleGroup.values):
        X_train, X_test = rpm1.iloc[train_index], rpm1.iloc[test_index]
        y_train, y_test = SampleGroup.values.ravel()[train_index], SampleGroup.values.ravel()[test_index]
        clf.fit(X_train, y_train)
        acc.append(accuracy_score(y_test, clf.predict(X_test)))
        auc.append(roc_auc_score(y_test, clf.predict(X_test)))
        sen.append(recall_score(y_test, clf.predict(X_test)))
        spe.append(precision_score(y_test, clf.predict(X_test)))
    print('acc:{}'.format(np.array(acc).mean()),
          'auc:{}'.format(np.array(auc).mean()),
          'spe:{}'.format(np.array(spe).mean()),
          'sen:{}'.format(np.array(sen).mean())
          )

# RF
def rf_cv():
    rfc = RFC()
    param_grid = {'n_estimators': [i for i in range(1, 101, 10)],
                  'max_depth': [i for i in range(10, 101, 10)],}
    grid = GridSearchCV(rfc, param_grid, cv=10, scoring='accuracy')
    grid.fit(rpm1, SampleGroup.values.ravel())
    print(grid.best_score_, grid.best_params_)

    rfc = RFC(n_estimators=11, max_depth=30)
    skf = StratifiedKFold(n_splits=10)
    acc = []
    auc = []
    spe = []
    sen = []
    for train_index, test_index in skf.split(rpm1, SampleGroup.values):
        X_train, X_test = rpm1.iloc[train_index], rpm1.iloc[test_index]
        y_train, y_test = SampleGroup.values.ravel()[train_index], SampleGroup.values.ravel()[test_index]
        rfc.fit(X_train, y_train)
        acc.append(accuracy_score(y_test, rfc.predict(X_test)))
        auc.append(roc_auc_score(y_test, rfc.predict(X_test)))
        sen.append(recall_score(y_test, rfc.predict(X_test)))
        spe.append(precision_score(y_test, rfc.predict(X_test)))
    print('acc:{}'.format(np.array(acc).mean()),
          'auc:{}'.format(np.array(auc).mean()),
          'spe:{}'.format(np.array(spe).mean()),
          'sen:{}'.format(np.array(sen).mean())
          )

# xgb
def xgb_cv():
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
    grid.fit(rpm1, SampleGroup.values.ravel())
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