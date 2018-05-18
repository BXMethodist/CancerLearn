import pandas as pd, numpy as np, os
from scipy.stats import norm
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from CancerROC import ROC_plot, ROC_plot2
from collections import defaultdict
from scipy import interp
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.metrics import accuracy_score

from sklearn import svm
from sklearn.multiclass import OneVsRestClassifier

def add_gene_length(target_df, gtf='../ref_data/hg19.ucscgenes.knowngene.xls', cell_type=None):
    gtf = pd.read_csv(gtf, sep='\t')
    gtf['hg19.kgXref.geneSymbol'] = gtf['hg19.kgXref.geneSymbol'].str.upper()

    length_results = []
    for i in target_df.index:
        cur_gtf = gtf[gtf['hg19.kgXref.geneSymbol'] == i]
        lengths = (cur_gtf['hg19.knownGene.txEnd'] - cur_gtf['hg19.knownGene.txStart']).max()
        length_results.append(lengths)

    target_df['gene_length'] = length_results

    return target_df

def get_columns(target1='tsg', target2='og', cutoff=-9):
    para_df = pd.read_excel('../grid/best_parameters.xlsx')
    para_df = para_df[(para_df['target1']==target1) & (para_df['target2']==target2)]
    if para_df.shape[0] < 1:
        para_df = para_df[(para_df['target1'] == target2) & (para_df['target2'] == target1)]
    para_df = para_df[para_df['logP']<cutoff]
    results = []
    for i in range(para_df.shape[0]):
        marker, feature, genebody = para_df.iloc[i, :3]
        results.append('_'.join([marker, feature, genebody]))
    return results

def label_label(table, label_table):
    result = [0]*table.shape[0]
    table['label'] = result
    print label_table['label'].unique()
    for i in label_table['label'].unique():
        cur_label = label_table[label_table['label']==i]
        table.ix[table.index.isin(cur_label.index), 'label'] = i
    return table

def center_normalization(table):
    return preprocessing.StandardScaler().fit(table)

def preprocessing_table(scaler, table):
    new_table_data = scaler.transform(table)
    new_table = pd.DataFrame(new_table_data, index=table.index, columns=table.columns)
    return new_table

def predict_LogisticRegression(X_train, Y_train, C=1., penalty='l1'):
    """
    :param table: a table which last columns is the label, (0 or 1), containing cell identity genes and non-cell identity genes
    :return: the sklearn logistic regression object
    """
    # print X_train
    predictor = LogisticRegression(penalty=penalty, C=C).fit(X_train, Y_train)
    return predictor

def predict_decision(predictor, real_table, label=False):
    if label:
        table = real_table.ix[:, :-1].copy()
    else:
        table = real_table.copy()
    result = predictor.decision_function(table)
    return pd.DataFrame(result, index=table.index)

def predict_proba(predictor, real_table, label=False):
    if label:
        table = real_table.ix[:, :-1].copy()
    else:
        table = real_table.copy()
    result = predictor.predict_proba(table)
    return pd.DataFrame(result, index=table.index)

def score(X, Y, predictor):
    return predictor.score(X, Y)

def ROC(y_true, y_score, title, labels):
    # # Compute ROC curve and ROC area for each class
    try:
        n_classes = y_true.shape[1]
    except:
        n_classes = 1
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[labels[i]], tpr[labels[i]], _ = roc_curve(y_true[:, i], y_score[:, i])

        roc_auc[labels[i]] = auc(fpr[labels[i]], tpr[labels[i]])
        # if pd.isnull(roc_auc[labels[i]]):
        #     print fpr[labels[i]], tpr[labels[i]], 'lala'
    return roc_auc, fpr, tpr

def get_training_testing(table, test_size=.2, random_state=0, n_class=2):
    Y = table[table.columns[-1]].tolist()
    # print table[table.columns[-1]]
    Y = label_binarize(Y, classes=range(n_class))
    # print Y
    X = table.iloc[:, :-1]
    # X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=random_state)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=random_state)
    # print X_test.shape
    return X_train, X_test, Y_train, Y_test

def all_AUC_attributes(CIG_df, train_df, groups, title, labels, test_size=0.2, C=.2, verbose=False, save=False, iterations=500, plot=False,
                       kind='train', method='AUC', band=False, fig=None, ax=None, color=None, test_function=get_training_testing,
                       train_cell_type_number=None, p='', penalty='l1', all_train_cell_type_number=None, hasLabel=False):
    # # #
    print 'title is ', title
    # train_df = train_df[['h3k4me3_total_width', 'h3k4me1_total_width', 'h3k27ac_total_width', 'h3k4me3_kurtosis']]
    # # #
    if not hasLabel:
        train_df = label_label(train_df, CIG_df)
    # # # #
    AUCs_train = {}
    TPRs_train = {}

    AUCs_test = {}
    TPRs_test = {}

    if method == 'score':
        scores = defaultdict(float)

    mean_fpr = np.linspace(0, 1, 101)

    original_df = None

    for r in range(iterations):
        # print r
        if all_train_cell_type_number is not None:
            if original_df is None:
                original_df = train_df.copy()
            cur_cell_types = original_df['cell_type'].sample(all_train_cell_type_number, random_state=r)
            train_df = original_df[original_df['cell_type'].isin(cur_cell_types)]

        if train_cell_type_number is None:
            X_train, X_test, Y_train, Y_test = test_function(train_df, random_state=r, test_size=test_size)
        else:
            X_train, X_test, Y_train, Y_test = test_function(train_df, random_state=r,
                                                             train_cell_type_number=train_cell_type_number)

        # if train_cell_type_number is None:
        #     X_train, X_test, Y_train, Y_test = test_function(train_df, random_state=r, test_size=test_size)
        # else:
        #     X_train, X_test, Y_train, Y_test = test_function(train_df, random_state=r,
        #                                                      train_cell_type_number=train_cell_type_number)
            # print len((Y_train[Y_train == 0]))
            # print (Y_train[Y_train == 1])
        # X_train, X_test, Y_train, Y_test = get_training_testing(train_df)
        # print Y_train.shape, X_train.shape

        y_trues_train = None
        y_scores_train = None

        y_trues_test = None
        y_scores_test = None

        first = True

        for i in range(len(groups)):
            cur_column = groups[i]
            cur_X_train = X_train.ix[:, cur_column]
            cur_X_test = X_test.ix[:, cur_column]

            if len(cur_X_train.shape) == 1:
                cur_X_train = cur_X_train.to_frame()
                cur_X_test = cur_X_test.to_frame()
            # print cur_X_train
            cur_CIG_predictor= predict_LogisticRegression(cur_X_train, Y_train, penalty=penalty, C=C)

            # print cur_CIG_predictor.coef_
            # print cur_X_train.columns
            if method == 'score':
                cur_score = score(cur_X_test, Y_test, cur_CIG_predictor)
                scores[labels[i]] += cur_score
            # print cur_column, score(cur_X_test, Y_test, cur_CIG_predictor)
            # print pd.DataFrame(cur_CIG_predictor.coef_.T, index=cur_X_train.columns)
            cur_y_score_train = predict_decision(cur_CIG_predictor, cur_X_train, False)
            cur_y_score_test = predict_decision(cur_CIG_predictor, cur_X_test, False)
            if first:
                y_trues_train = Y_train
                y_trues_test = Y_test
                y_scores_train = cur_y_score_train.values
                y_scores_test = cur_y_score_test.values
                first = False
            else:
                y_trues_train = np.concatenate((y_trues_train, Y_train), axis=1)
                y_trues_test = np.concatenate((y_trues_test, Y_test), axis=1)
                y_scores_train = np.concatenate((y_scores_train, cur_y_score_train.values), axis=1)
                y_scores_test = np.concatenate((y_scores_test, cur_y_score_test.values), axis=1)

        # pd.DataFrame(np.concatenate([y_trues_test, y_scores_test], axis=1)).to_csv(title+'_for_prism.csv', index=None, columns=None)

        auc_train, fpr_train, tpr_train = ROC(y_trues_train, y_scores_train, 'train_'+title, labels)
        auc_test, fpr_test, tpr_test = ROC(y_trues_test, y_scores_test, 'test_'+title, labels)

        for label in labels:
            if kind == 'train' and band:
                plt.plot(fpr_train[label], tpr_train[label], lw=0.4, alpha=0.1, color='grey')
            elif kind == 'test' and band:
                plt.plot(fpr_test[label], tpr_test[label], lw=0.4, alpha=0.1, color='grey')
            tpr_train[label] = interp(mean_fpr, fpr_train[label], tpr_train[label])
            tpr_test[label] = interp(mean_fpr, fpr_test[label], tpr_test[label])


        for l in range(len(labels)):
            label = labels[l]
            if label not in AUCs_train:
                AUCs_train[label] = auc_train[label]
                AUCs_test[label] = auc_test[label]
                TPRs_train[label] = [tpr_train[label]]
                TPRs_test[label] = [tpr_test[label]]
            else:
                AUCs_train[label] += auc_train[label]
                AUCs_test[label] += auc_test[label]
                TPRs_train[label].append(tpr_train[label])
                TPRs_test[label].append(tpr_test[label])
    for label in labels:
        AUCs_train[label] /= iterations
        AUCs_test[label] /= iterations

        if method=='score':
            scores[label] /= iterations

    if plot:
        if kind == 'train':
            ROC_plot2(AUCs_train, mean_fpr, TPRs_train, len(labels), title, labels, verbose=verbose, save=save, band=band, fig=fig, ax=ax, color=color, p=p)
        elif kind == 'test':
            ROC_plot2(AUCs_test, mean_fpr, TPRs_test, len(labels), title, labels, verbose=verbose, save=save, band=band, fig=fig, ax=ax, color=color, p=p)

    if method == 'AUC':
        results_train_df = pd.DataFrame.from_dict(AUCs_train, orient='index')
        results_train_df.columns = ['train']
        results_test_df = pd.DataFrame.from_dict(AUCs_test, orient='index')
        results_test_df.columns = ['test']

        result_df = results_train_df.join(results_test_df)
        return result_df
    elif method == 'score':
        result_df = pd.DataFrame.from_dict(scores, orient='index')
        result_df.columns = ['accuracy']
    elif method == 'cutoff':
        cutoffs = {}
        for label in labels:
            cutoffs[label] = np.sqrt(np.min((0 - mean_fpr)**2 + (1.0-np.asarray(np.mean(TPRs_test[label], axis=0)))**2))

            # plt.plot(mean_fpr, np.mean(TPRs_test[label], axis=0))
            # plt.show()
            # print cutoffs[label], label
        result_df = pd.DataFrame.from_dict(cutoffs, orient='index')
        result_df.columns = ['distance']

    if kind =='test':
        for key in TPRs_test.keys():
            TPRs_test[key] = np.mean(TPRs_test[key], axis=0)

        TPRs_test_df = pd.DataFrame.from_dict(TPRs_test)
        TPRs_test_df.to_csv(title+'TPR.csv')

    return result_df



if __name__ == '__main__':

    ## predict the lncRNA with CIG model and get the top list
    if True:
        # og_df1 = pd.read_csv('../genelist/OG_top500.xls', sep='\t', index_col=0)
        # og_df1['label'] = [0]*og_df1.shape[0]
        #
        # tsg_df1 = pd.read_csv('../genelist/TSG_top500.xls', sep='\t', index_col=0)
        # tsg_df1['label'] = [1] * tsg_df1.shape[0]
        #
        # og_genes = set(og_df1.index)
        # tsg_genes = set(tsg_df1.index)
        # #
        # cancer = []
        #
        # print len(og_genes), len(tsg_genes)
        #
        # for o in og_genes:
        #     if o not in tsg_genes:
        #         cancer.append([o, 0])
        # for t in tsg_genes:
        #     if t not in og_genes:
        #         cancer.append([t, 1])
        # # print og_df.shape, tsg_df.shape
        #
        # cancer_df = pd.DataFrame(cancer)
        # cancer_df.columns=['gene', 'label']
        # cancer_df = cancer_df.set_index(['gene'])
        # cancer_df['cell_type'] = ['CD4']*cancer_df.shape[0]

        cancer_df = pd.read_excel('../genelist/TUSON_OGTSG_subgroups.xlsx', index_col=0)
        cancer_df = cancer_df.iloc[:, [0,1]]
        cancer_df.columns = ['type', 'label']

        train_df = pd.read_csv('../train/Cancer_CD4_default.csv', index_col=0)
        train_df.index = train_df.index.str.upper()

        train_df = train_df[train_df.index.isin(cancer_df.index)]

        # columns = get_columns(target1='tsg', target2='og', cutoff=-10)
        columns = [x for x in train_df.columns if x.find('single') == -1]

        train_df = train_df[columns]

        train_df = train_df.dropna()
        train_df = train_df.drop_duplicates()

        train_df = add_gene_length(train_df)

        """
        add information from other cell types
        """
        if True:
            tables = ['../identity/' + x for x in os.listdir('../identity/') if
                      x.find('ESC') == -1 and x.endswith('.csv')]
            cell_types = [x.split('_')[0] for x in os.listdir('../identity') if
                          x.find('ESC') == -1 and x.endswith('.csv')]
            columns = pd.read_csv(tables[0], index_col=0).columns
            columns = [x for x in columns if
                       x.find('single') == -1 and x.find('length') == -1 and x.find('label') == -1]

            for c in columns:
                other_og_df = pd.read_excel('../identity/'+'OG500_'+c+'.xlsx', index_col=0)
                other_tsg_df = pd.read_excel('../identity/' + 'TSG500_' + c + '.xlsx', index_col=0)
                print other_og_df.shape, other_tsg_df.shape
                other_df = other_og_df.append(other_tsg_df)
                other_df.columns =[x+'_'+c for x in other_df.columns]
                print c
                print train_df.shape, other_df.shape

                train_df = pd.concat([train_df, other_df], axis=1)



        train_df = label_label(train_df, cancer_df)



        train_df.to_csv('../train/Cancer_train_label.xls', sep='\t')






        scaler = center_normalization(train_df.iloc[:, :-1])

        train_df.iloc[:, :-1] = preprocessing_table(scaler, train_df.iloc[:, :-1])


        model = LogisticRegression(penalty='l1', C=.1)
        model = OneVsRestClassifier(svm.SVC(probability=True, C=.2))
        # ROC_plot(model, train_df.iloc[:, :-1], train_df.iloc[:, -1])
        ROC_plot(model, train_df.iloc[:, :-1], train_df.iloc[:, -1], classes=[0,1,2,3])

        # model = OneVsRestClassifier(svm.SVC(probability=True))
        # model.fit(train_df.iloc[:, :-1], train_df.iloc[:, -1])
        # result = predict_decision(model, train_df.iloc[:, :-1])
        # result2 = predict_proba(model, train_df.iloc[:, :-1])
        #
        # result = np.concatenate([result, result2], axis=1)
        # result_df = pd.DataFrame(result, index=train_df.index)
        # result_df['predict_label'] = model.predict(train_df.iloc[:, :-1])
        #
        # result_df['label'] = train_df.ix[result_df.index, 'label']

        # result_df = result_df[result_df['label']==1]
        # result_df = result_df[result_df['label'] != result_df['predict_label']]
        #
        # result_df.to_csv('../results/training_prediction.csv')


