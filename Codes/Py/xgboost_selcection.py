import pandas as pd
import re
import numpy as np
import math
import xgboost as xgb
from xgboost import XGBClassifier
from xgboost import plot_importance
from matplotlib import pyplot
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel

def readFiles(filepath_counts,filepath_names):
    counts_data = pd.read_csv(filepath_counts, sep = ',', index_col = 0)
    DEGnames =  list(pd.read_csv(filepath_names, sep = ','
                                 ,index_col = False).iloc[:,0])
    
    return counts_data, DEGnames

def dataProcessing(counts_data, DEGnames):
    RNA_names = counts_data.index.tolist()
    sample_name = counts_data.columns.tolist()
    sample_name_l = []
    
    for gene in DEGnames:
        if gene in RNA_names:
            continue
        else:
            counts_data = counts_data.drop(gene)
    
    for i in sample_name:
        if re.match(r'^TCGA\.(.*)\.(.*)\.01A\.(.*)\.(.*)',i):
            i = 1
            sample_name_l.append(i)
        else:
            i = 0
            sample_name_l.append(i)

    counts_data.columns = sample_name_l
    counts_data = counts_data.T
    counts_data.index.name = 'Type'
    counts_data.replace(0, 1,inplace = True)
   
    f = lambda num:math.log(num + 1, 10)
    counts_data = counts_data.applymap(f).reset_index()
    
    return counts_data


def fea_output(counts_data, DEGnames, times):
    X=counts_data.iloc[:,1:].values
    y=counts_data.iloc[:,0].values
    
    #feat_labels=pd.DataFrame(miRNA_DEG.columns[1:])
    result_importance = pd.DataFrame(index = DEGnames)
    for i in range(times):
        x_train, x_test, y_train, y_test = train_test_split(X, y
                                                        , random_state = i
                                                        , stratify=y
                                                        , test_size=0.3)
                                                   
                                                      
        #feature_names = list(miRNA_DEG.drop(['Type'], axis=1).columns)	
        # 得到所有的特征

        model = xgb.XGBClassifier(eta=0.01, gamma=1, learning_rate=0.1
                                  , max_depth=3
                                  , n_jobs=-1
                                  , importance_type='total_gain').fit(x_train,
                                                                      y_train)

        #整合特征100次循环的重要性打分
        importances1 = model.feature_importances_    
        result_importance['run_%d' % i] = importances1
       # importances = list(model.feature_importances_)
        #importances_df = pd.DataFrame(importances)
    thred = np.mean(result_importance)
    result_importance1 = np.sum(result_importance > thred, axis=1)
    fea_sel = feat_labels.loc[result_importance1.values > 50,:].values
    fea_sel_score_df = result_importance.loc[fea_sel[:,0], :]
    fea_sel_score_df.to_csv(r"selection_xgboost.csv"
                            , header = True
                            , index = True)

    
    
    