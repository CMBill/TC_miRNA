# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 10:18:26 2021

@author: 12870
"""
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.metrics import f1_score, matthews_corrcoef
from sklearn.metrics import accuracy_score, precision_score, recall_score
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, plot_roc_curve, auc
import matplotlib.pyplot as plt

data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\generank_selected.csv')
X, y = np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0])

X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,
                                                    test_size=0.3,
                                                    random_state=5)
model_ = svm.SVC(kernel='linear', C=1)
model_.fit(X_train, y_train)
print(metrics.roc_auc_score(y_test, model_.predict(X_test)))

import numpy as np
import pandas as pd
from sklearn import svm
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score, matthews_corrcoef
from sklearn.linear_model import LogisticRegression
import xgboost as xgb

if __name__ == "__main__":
    path = 'D:\\project\\TCGA_stomach_analysis\\data\\generank_selected.csv'  # 数据文件路径
    data = pd.read_csv(path, index_col=False)
    x = np.array(data.iloc[:, 1:], dtype=np.float)
    y = np.array(data.iloc[:, 0], dtype=np.float)
    # 即x为miRNA里的前41列特征列，y为miRNA里的第42列标签列
    x_train, x_test, y_train, y_test = train_test_split(x, y, random_state=100, stratify=y, train_size=0.7,
                                                        test_size=0.3)  #

    # 分类器
    # clf = svm.SVC(C=0.1, kernel='rbf', decision_function_shape='ovr')
    # ovr   one vs rest
    clf = svm.SVC(kernel='linear')
    rfc = RFC(n_estimators=75, max_depth=12, min_samples_split=2).fit(X_train,
                                                                      y_train)
    clf.fit(x_train, y_train.ravel())
    lr = LogisticRegression(penalty='l2', random_state=0, C=1).fit(X_train, y_train)
    XGB = xgb.XGBClassifier().fit(X_train, y_train)
print(clf.score(x_train, y_train))  # 精度
print('训练集准确率：', accuracy_score(y_train, clf.predict(x_train)))
print('训练集精度：', precision_score(y_train, clf.predict(x_train)))
print(clf.score(x_test, y_test))
print('测试集acc：', accuracy_score(y_test, clf.predict(x_test)))
print('测试集SPE：', precision_score(y_test, clf.predict(x_test)))
print('测试集AUC：', metrics.roc_auc_score(y_test, clf.predict(x_test)))
print('测试集f1:', f1_score(y_test, clf.predict(x_test)))
print('测试集sen:', recall_score(y_test, clf.predict(x_test)))

print('测试集acc：', accuracy_score(y_test, rfc.predict(x_test)))
print('测试集SPE：', precision_score(y_test, rfc.predict(x_test)))
print('测试集AUC：', metrics.roc_auc_score(y_test, rfc.predict(x_test)))
print('测试集f1:', f1_score(y_test, rfc.predict(x_test)))
print('测试集sen:', recall_score(y_test, rfc.predict(x_test)))

print('测试集acc：', accuracy_score(y_test, lr.predict(x_test)))
print('测试集SPE：', precision_score(y_test, lr.predict(x_test)))
print('测试集AUC：', metrics.roc_auc_score(y_test, lr.predict(x_test)))
print('测试集f1:', f1_score(y_test, lr.predict(x_test)))
print('测试集sen:', recall_score(y_test, lr.predict(x_test)))

print('测试集acc：', accuracy_score(y_test, XGB.predict(x_test)))
print('测试集SPE：', precision_score(y_test, XGB.predict(x_test)))
print('测试集AUC：', metrics.roc_auc_score(y_test, XGB.predict(x_test)))
print('测试集f1:', f1_score(y_test, XGB.predict(x_test)))
print('测试集sen:', recall_score(y_test, XGB.predict(x_test)))
