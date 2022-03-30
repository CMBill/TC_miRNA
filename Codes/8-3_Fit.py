import pandas as pd
import numpy as np
from sklearn import svm
import xgboost as XGB
from sklearn.ensemble import RandomForestClassifier as RFC
import Codes.fit_fc as fc

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

print('svm')
clf = svm.SVC(kernel='linear', C=1, probability=True)
fc.fit_c(clf, X_rfe, y_rfe, 5)
print('rf')
rfc = RFC(n_estimators=61, max_depth=70)
fc.fit_c(rfc, X_rfe, y_rfe, 5)
print('xgb')
xgb = XGB.XGBClassifier(objective='binary:logistic',
                        nthread=4,
                        learning_rate=0.1,
                        n_estimators=170,
                        max_depth=3,
                        min_child_weight=1,
                        gamma=0.0,
                        colsample_bytree=0.6,
                        subsample=0.6,
                        reg_alpha=0.01,
                        reg_lambda=1)
fc.fit_c(xgb, X_rfe, y_rfe, 5)
