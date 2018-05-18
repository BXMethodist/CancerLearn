print(__doc__)

import matplotlib.pyplot as plt, numpy as np, pandas as pd
from itertools import cycle
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from scipy import interp
from sklearn import preprocessing
from sklearn.metrics import accuracy_score

def ROC_plot(model, X, y, fname, classes=[0,1,2], test_size=.1, iterations=100, verbose=False):
    # Import some data to play with
    # Binarize the output
    y = label_binarize(y, classes=classes)
    n_classes = y.shape[1]
    # n_classes = len(classes)
    # shuffle and split training and test sets
    # Learn to predict each class against the other
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    mean_fpr = np.linspace(0, 1, 101)
    for r in range(iterations):
        # try:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=r)
        predictor = model.fit(X_train, y_train)

        # y_score = predictor.decision_function(X_test)
        y_score = predictor.predict_proba(X_test)
        for i in range(n_classes):
            cur_fpr, cur_tpr, _ = roc_curve(y_test[:, i], y_score[:, i])
            if pd.isnull(np.sum(cur_tpr)):
                print cur_tpr
                continue

            cur_tpr = interp(mean_fpr, cur_fpr, cur_tpr)

            if i in tpr.keys():
                tpr[i] += cur_tpr
                roc_auc[i] += auc(mean_fpr, cur_tpr)
            else:
                tpr[i] = cur_tpr
                roc_auc[i] = auc(mean_fpr, cur_tpr)
        # except:
        #     print 'fail'
        #     continue
    # print tpr
    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(mean_fpr)
    for i in range(n_classes):
        mean_tpr += interp(mean_fpr, mean_fpr, tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= n_classes

    fpr["macro"] = mean_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
    labels = ['Control', 'TSG', 'OG']
    colors = cycle(['black', 'darkorange', 'green', 'aqua', 'cornflowerblue', 'purple'])
    for i, color in zip(range(n_classes), colors):
        if i == 0:
            continue
        plt.plot(mean_fpr, tpr[i]/iterations, color=color,
                 label='ROC curve of group {0} (area = {1:0.2f})'
                 ''.format(str(labels[i]), roc_auc[i]/iterations))
        pass
    plt.plot(mean_fpr, tpr['macro'] / iterations, color='red',
             label='ROC curve of {0} (area = {1:0.2f})'
                   ''.format('full model', roc_auc['macro']/iterations))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig('../plots/'+fname, transparent=True)



def ROC_plot2(roc_auc, fpr, tpr, n_classes, title, labels, verbose=False, save=False, band=False, plot=False, fig=None, ax=None, color=None,p='', legend=False):
    # plt.figure()
    if not plot:
        if fig is None and ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    color_idx = np.linspace(0, 1, n_classes)
    colors = ['red','black', 'cyan', 'pink', 'orange', 'skyblue', 'yellow', 'pink', 'purple',]
    # colors = ['black','cyan', 'red', 'black', 'grey', 'green', 'skyblue', 'gray', 'cyan', 'purple',]
    # colors = ['red', 'blue', 'grey', 'grey', 'green', 'skyblue', 'gray', 'cyan', 'purple', ]
    c=0
    font = {'fontname': 'Helvetica', 'fontsize': 15}
    for i, j in zip(range(n_classes), color_idx):
        print i, j
        y_info = np.mean(tpr[labels[i]], axis=0)
        y_info[0] =0
        if color is None:
            plt.plot(fpr[0:], y_info, color=colors[c], linestyle='-', alpha=1,
                     label='{0} (area = {1:0.2f})'
                           ''.format(labels[i], roc_auc[labels[i]]))
            c += 1
        else:
            plt.plot(fpr[0:], y_info[0:], color=color,linestyle='--',alpha=1,
                     label='{0} (area = {1:0.2f})'
                           ''.format(labels[i], roc_auc[labels[i]]))
            # label = '{0} (area = {1:0.2f})'
            # ''.format(labels[i], roc_auc[labels[i]])

        if band:
            mean_tpr = np.mean(tpr[labels[i]], axis=0)
            mean_tpr[-1] = 1.0
            std_tpr = np.std(tpr[labels[i]], axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            plt.fill_between(fpr, tprs_lower, tprs_upper, color=plt.cm.plasma(j), alpha=.3)
    if p != '':
        ax.text(0.95, 0.5, 'p value: '+p, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes,
                **font)

    # plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', **font)
    plt.ylabel('True Positive Rate', **font)
    plt.title(title)
    if legend:
        plt.legend(loc="lower right")
    if verbose:
        plt.show()
    # plt.axis('off')
    if save:
        if title == '':
            plt.savefig('result' + '.pdf', transparent=True)
        else:
            plt.savefig(title+'.pdf', transparent=True)

        plt.close('all')
    return roc_auc, fpr, tpr
