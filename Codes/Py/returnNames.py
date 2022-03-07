# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 16:59:12 2021

@author: LanHao
"""


def readNames(i):
    import os
    os.chdir('D:/project/TCGA_stomach_analysis/data')
    filepath_counts = ['./download_mirna/miRNA_expr1.csv',
                       './download_rnaseq/mRNA/mRNA_expr.csv',
                       './download_rnaseq/lncRNA/lncRNA_expr.csv']

    filepath_names = ['./download_mirna/miRNA_DEGname.csv',
                      './download_rnaseq/mRNA/mRNA_DEGname.csv',
                      './download_rnaseq/lncRNA/lncRNA_DEGname.csv']
    return filepath_counts[i], filepath_names[i]


def saveNames(i, j):
    import os
    os.chdir('D:/project/TCGA_stomach_analysis/data')
    name_ls = ['miRNA', 'mRNA', 'lncRNA']
    miRNA = ['./download_mirna/miRNA_lr.csv',
             './download_mirna/miRNA_knn.csv',
             './download_mirna/miRNA_svm.csv',
             './download_mirna/miRNA_rf.csv',
             './download_mirna/miRNA_xgboost.csv',
             './download_mirna/miRNA_dt.csv']
    mRNA = ['./download_rnaseq/mRNA/mRNA_lr.csv',
            './download_rnaseq/mRNA/mRNA_knn.csv',
            './download_rnaseq/mRNA/mRNA_svm.csv',
            './download_rnaseq/mRNA/mRNA_rf.csv',
            './download_rnaseq/mRNA/mRNA_xgboost.csv',
            './download_rnaseq/mRNA/mRNA_dt.csv']
    lncRNA = ['./download_rnaseq/lncRNA/lncRNA_lr.csv',
              './download_rnaseq/lncRNA/lncRNA_knn.csv',
              './download_rnaseq/lncRNA/lncRNA_svm.csv',
              './download_rnaseq/lncRNA/lncRNA_rf.csv',
              './download_rnaseq/lncRNA/lncRNA_xgboost.csv',
              './download_rnaseq/lncRNA/lncRNA_dt.csv']
    filepath_save = dict(miRNA=miRNA, mRNA=mRNA, lncRNA=lncRNA)

    return filepath_save[name_ls[i]][j]
