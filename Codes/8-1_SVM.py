import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn import svm
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score, confusion_matrix, roc_auc_score
from sklearn.metrics import recall_score
import xgboost as xgb

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
SampleGroup.iloc[SampleGroup.iloc[:, 1] == 'cancer', 1] = 1
SampleGroup.iloc[SampleGroup.iloc[:, 1] == 'normal', 1] = 0

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
x_train, x_test, y_train, y_test = train_test_split(rpm1, SampleGroup, random_state=100, stratify=SampleGroup,
                                                    train_size=0.7, test_size=0.3)
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
    f1 = []
    spe = []
    sen = []
    for train_index, test_index in skf.split(rpm1, SampleGroup.values.reshape(-1,1)):
        X_train, X_test = rpm1.iloc[train_index], rpm1.iloc[test_index]
        y_train, y_test = SampleGroup.values.ravel()[train_index], SampleGroup.values.ravel()[test_index]
        clf.fit(X_train, y_train)
        cm = confusion_matrix(y_true=y_test, y_pred=clf.predict(X_test))
        print(cm)
        tn, fp, fn, tp = cm.ravel()
        acc.append((tp + tn) / (tp + tn + fp + fn))
        auc.append(roc_auc_score(y_test, clf.predict(X_test)))
        f1.append(f1_score(y_test, clf.predict(X_test)))
        sen.append(tp / (tp + fn))
        spe.append(tn / (tn + fp))
    print('acc:{}'.format(np.array(acc).mean()),
          'auc:{}'.format(np.array(auc).mean()),
          'f1:{}'.format(np.array(f1).mean()),
          'spe:{}'.format(np.array(spe).mean()),
          'sen:{}'.format(np.array(sen).mean())
          )
