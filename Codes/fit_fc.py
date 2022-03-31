import numpy as np
from matplotlib import pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, roc_curve, auc, RocCurveDisplay
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score, recall_score, f1_score
import time


def metrics_fit(classifier, X, y, k):
    """
    """
    skf = StratifiedKFold(n_splits=k)
    acc = []
    auc1 = []
    f1 = []
    spe = []
    sen = []
    tprs = []
    aucs = []
    plt.figure(figsize=(10, 10))
    mean_fpr = np.linspace(0, 1, 100)
    i = 0
    for train, test in skf.split(X, y):
        classifier.fit(X[train], y[train].ravel())
        probas_ = classifier.predict_proba(X[test])
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1], pos_label='1')
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.3,
                 label='ROC fold %d (AUC = %0.4f)' % (i, roc_auc))
        i += 1
        acc.append(accuracy_score(y[test], classifier.predict(X[test])))
        auc1.append(roc_auc_score(y[test], classifier.predict(X[test])))
        f1.append(f1_score(y[test], classifier.predict(X[test]), pos_label='1'))
        sen.append(recall_score(y[test], classifier.predict(X[test]), pos_label='1'))
        spe.append(precision_score(y[test], classifier.predict(X[test]), pos_label='1'))
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Chance', alpha=.8)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.4f )' % mean_auc,
             lw=2, alpha=.8)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.xlabel('False Positive Rate', fontsize=18)
    plt.ylabel('True Positive Rate', fontsize=18)
    plt.title('Cross-Validation ROC of SVM', fontsize=18)
    plt.legend(loc="lower right", prop={'size': 15})
    print('acc:{}'.format(np.array(acc).mean()),
          'auc:{}'.format(np.array(auc1).mean()),
          'f1:{}'.format(np.array(f1).mean()),
          'spe:{}'.format(np.array(spe).mean()),
          'sen:{}'.format(np.array(sen).mean())
          )


def fit_c(classifier, X, y, k):
    """
    :param classifier: 模型及参数
    :param X: 表达值
    :param y: 标签
    :param k: K折数
    :return: 输出ROC曲线及评价参数，无返回值
    """
    cv = StratifiedKFold(n_splits=k)
    acc = []
    auc1 = []
    f1 = []
    spe = []
    sen = []
    tprs = []
    aucs = []
    # plt.figure(figsize=(10, 10))
    plt.rcParams['figure.figsize'] = (10.0, 6.0)
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(cv.split(X, y)):
        classifier.fit(X[train], y[train].ravel())
        viz = RocCurveDisplay.from_estimator(classifier, X[test], y[test],
                                             name='ROC fold {}'.format(i),
                                             alpha=0.3, lw=1, ax=ax)
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
        acc.append(accuracy_score(y[test], classifier.predict(X[test])))
        auc1.append(roc_auc_score(y[test], classifier.predict(X[test])))
        f1.append(f1_score(y[test], classifier.predict(X[test]), pos_label='1'))
        sen.append(recall_score(y[test], classifier.predict(X[test]), pos_label='1'))
        spe.append(precision_score(y[test], classifier.predict(X[test]), pos_label='1'))

    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.4f $\pm$ %0.4f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2)

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
           title="Receiver operating characteristic")
    ax.legend(loc="lower right")
    t = time.gmtime()
    timestamp = time.strftime('%m%d%H%M%S', t)
    plt.savefig('./pics/Model' + timestamp + '.pdf')
    plt.show()

    print('./pics/Model' + timestamp + '.pdf' + '\n' +
          'acc:{}'.format(np.array(acc).mean()),
          'auc:{}'.format(np.array(auc1).mean()),
          'f1:{}'.format(np.array(f1).mean()),
          'spe:{}'.format(np.array(spe).mean()),
          'sen:{}'.format(np.array(sen).mean())
          )
