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

print('svm:')
clf = svm.SVC(kernel='linear', C=1, probability=True)  # kernel='linear', C=1
fc.fit_c(clf, X_rfe, y_rfe, 5)
print('rf:')
rfc = RFC(n_estimators=61, max_depth=20)
fc.fit_c(rfc, X_rfe, y_rfe, 5)
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
fc.fit_c(xgb, X_rfe, y_rfe, 5)
