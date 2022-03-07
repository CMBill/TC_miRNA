# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 09:09:28 2021

@author: 12870
"""

from sklearn.metrics import roc_curve, plot_roc_curve, auc
from model_evaluation import auc_, mcc_, f1_, precision_, sen_, acc_
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn import svm
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE
import numpy as np
import pandas as pd
from sklearn import metrics
data = pd.read_csv('D:\\project\\TCGA_stomach_analysis\\data\\generank_selected.csv')
X, y = np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0])

X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    test_size=0.3,
                                                    random_state=0)

model_ = svm.SVC(probability=True)
model_.fit(X_train, y_train)
prob = model_.predict(X_test)
fpr, tpr, thresholds = roc_curve(y_test, prob)
roc_auc = auc(fpr, tpr)
print(metrics.roc_auc_score(y_test,model_.predict(X_test)))
def drawRoc(roc_auc,fpr,tpr):
    plt.subplots(figsize=(7, 5.5))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.show()
drawRoc(roc_auc, fpr, tpr)
