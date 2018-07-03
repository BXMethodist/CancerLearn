"""
This module is designed to take the environment context to change the estimated probablity for OG and TSG
"""

import pandas as pd, numpy as np, os
from scipy import stats
from sklearn import preprocessing
from sklearn.preprocessing import label_binarize
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MinMaxScaler

eps = np.finfo(float).eps

def sigmoid(x):
    return (1 / (1 + np.exp(-x)))

def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

def label_label(table, positive_table):
    result = [0]*table.shape[0]
    table['label'] = result
    table.loc[table.index.isin(positive_table.index), 'label'] = 1
    return table

def center_normalization(table):
    return preprocessing.StandardScaler().fit(table)

def preprocessing_table(scaler, table):
    new_table_data = scaler.transform(table)
    new_table = pd.DataFrame(new_table_data, index=table.index, columns=table.columns)
    return new_table

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

def predict_LogisticRegression(X_train, Y_train, C=1., penalty='l1'):
    """
    :param table: a table which last columns is the label, (0 or 1), containing cell identity genes and non-cell identity genes
    :return: the sklearn logistic regression object
    """
    # print X_train
    predictor = LogisticRegression(penalty=penalty, C=C).fit(X_train, Y_train)
    return predictor


if __name__ == "__main__":
    """
    predict the context with default parameters.
    """
    #### predict HUVEC get the target genes in top rank.
    if False:
        CIG_df = pd.read_excel('../ref_data/CIG_Bo_curated_311.xlsx', index_col=0)
        # CIG_df = pd.read_csv('CIG_strong_bo.csv', index_col=0)
        results = []
        # for r in range(0, 10000, 100):
        for r in [100]:
            train_df = pd.read_csv('../context/Cancer_Learn_CIG_train.csv', index_col=0)
            columns = [c for c in train_df.columns if (c.find('27ac')!=-1 or c.find('4me1')!=-1 or c.find('4me3')!=-1) and c.find('cell_type')==-1 and (c.find('total_width_genebody')!=-1 or c.find('height_genebody')!=-1)]
            train_df = train_df[columns]

            train_df = label_label(train_df, CIG_df)
            #
            scaler = center_normalization(train_df.iloc[:, :-1])
            # print scaler.mean_.shape
            #
            train_df.iloc[:, :-1] = preprocessing_table(scaler, train_df.iloc[:, :-1])

            predictor = predict_LogisticRegression(train_df.iloc[:, :-1], train_df.iloc[:, -1], C=.2)

            cell_types = ['../context/'+x for x in os.listdir('../context') if x.endswith('real_table.csv')]
            # cell_types = ['HUVEC']

            from scipy.stats import norm
            from rpy2.robjects.packages import importr
            from rpy2.robjects.vectors import FloatVector

            for c in cell_types:
                celltype = c.replace('../context/','').split('_')[0]
                real_table = pd.read_csv(c, dtype={'gene_id': str})
                real_table = real_table.set_index('gene_id')

                real_table = real_table[columns]
                real_table = preprocessing_table(scaler, real_table)

                result = predict_decision(predictor, real_table)
                result2 = predict_proba(predictor, real_table)
                # print result.shape, result2.shape
                #
                result = np.concatenate([result, result2], axis=1)
                #
                result_df = pd.DataFrame(result, index=real_table.index)
                result_df.columns = ['distance', 'non-CIG_prob', 'CIG_prob']
                del result_df['non-CIG_prob']
                # del result_df['CIG_prob']
                result_df.to_csv('../context/'+celltype + '_all_features_result.csv')

                stats = importr('stats')
                df = pd.read_csv('../context/'+ celltype + '_all_features_result.csv')

                df["p_value"] = 1 - norm.cdf(df['distance'], scale=1.25)
                df['FDR'] = stats.p_adjust(FloatVector(df["p_value"].tolist()),
                                           method='BH')

                # print df[df['p_value'] < 0.01].shape[0]
                df = df.sort_values(by=['distance'], ascending=False)
                # df['rank'] = np.arange(0, 1, 1./df.shape[0])
                # df['rank'] = 1./(1+(0.05/(df['rank']+eps)))**4
                # df['rank'] = 1-MinMaxScaler().fit(df['rank'].values.reshape([-1,1])).transform(df['rank'].values.reshape([-1,1]))

                # del df['CIG_prob']
                # del df['non-CIG_prob']
                df = df.set_index(['gene_id'])

                df.to_csv('../context/'+ celltype + '_all_features_result.csv')


    """
    add environment context to do the prediction
    """
    if False:
        # known_gtf = pd.read_csv('../ref_data/hg19.ucscgenes.knowngene.xls', sep='\t')
        known_gtf = pd.read_csv('../ref_data/hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')

        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)

        pan_og_candidates = cancer_df[cancer_df['type']=='OG'].index
        pan_tsg_candidates = cancer_df[cancer_df['type']=='TSG'].index

        predict_df = pd.read_excel('../results/aggregate_5_vote100.xlsx', index_col=0)
        predict_df = predict_df[~predict_df.index.duplicated(keep='first')]

        exp_df = pd.read_csv('../data/breast/CCLE_breast.xls', index_col=0, sep='\t')
        exp_df = exp_df[~exp_df.index.duplicated(keep='first')]
        exp_df = quantileNormalize(exp_df)

        exp_columns = ['HMEL_BREAST', 'MCF7_BREAST', 'MDAMB231_BREAST', 'AU565_BREAST', 'ZR751_BREAST',
                       'MDAMB361_BREAST', 'UACC812_BREAST', 'SKBR3_BREAST',
                       'HCC1954_BREAST', 'MDAMB468_BREAST', 'HCC1937_BREAST', 'MDAMB436_BREAST']
        celltype_columns = ['HMEL', 'MCF7', 'MB231', 'AU565', 'ZR751',
                       'MB361', 'UACC812', 'SKBR3',
                       'HCC1954', 'MB468', 'HCC1937', 'MB436']

        # exp_columns = ['HMEL_BREAST', 'MCF7_BREAST', 'MDAMB231_BREAST']
        # celltype_columns = ['HMEL', 'MCF7', 'MB231']

        result_df = pd.DataFrame(index=celltype_columns, columns=['p_og', 'p_tsg', 'pan_p_og', 'pan_p_tsg'])
        result_df.index.name = 'cell_type'

        for i in range(1, len(exp_columns)):
            cur_exp_column = exp_columns[i]
            cur_celltype_column = celltype_columns[i]

            cur_context_df = pd.read_csv('../context/'+cur_celltype_column+'_all_features_result.csv', index_col=0)
            cur_context_df = cur_context_df[~cur_context_df.index.duplicated(keep='first')]

            common_index = set(predict_df.index).intersection(cur_context_df.index).intersection(exp_df.index)

            # print len(common_index)

            cur_predict_df = predict_df[predict_df.index.isin(common_index)].copy()
            cur_context_df = cur_context_df[cur_context_df.index.isin(common_index)]
            exp_df = exp_df[exp_df.index.isin(common_index)]

            cur_distance = cur_context_df['distance'].values.reshape(-1, 1)
            cur_context_df['og_distance'] = MinMaxScaler().fit(cur_distance).transform(cur_distance)
            cur_context_df['tsg_distance'] = MinMaxScaler().fit(cur_distance*-1).transform(cur_distance*-1)

            # cur_predict_df[cur_celltype_column + '_OG_prob'] = cur_predict_df['OG_prob'] + cur_context_df['distance']
            cur_predict_df[cur_celltype_column + '_OG_prob'] = cur_predict_df['OG_prob'] * cur_context_df['og_distance'] + cur_context_df['og_distance']
            cur_predict_df[cur_celltype_column + '_TSG_prob'] = (cur_predict_df['TSG_prob']+cur_predict_df['OG_prob']) * cur_context_df['tsg_distance'] - cur_predict_df['OG_prob']**10 + cur_context_df['tsg_distance']

            cur_predict_df = cur_predict_df.sort_values(by=[cur_celltype_column + '_OG_prob'], ascending=False)
            cur_predict_df[cur_celltype_column + '_OG_prob'] = np.arange(0, 1, 1./cur_predict_df.shape[0])
            cur_predict_df[cur_celltype_column + '_OG_prob'] = 1. / (1 + (0.05 / (cur_predict_df[cur_celltype_column + '_OG_prob'] + eps))) ** 4
            cur_predict_df[cur_celltype_column + '_OG_prob'] = 1 - MinMaxScaler().fit(cur_predict_df[cur_celltype_column + '_OG_prob'].values.reshape([-1, 1])).transform(cur_predict_df[cur_celltype_column + '_OG_prob'].values.reshape([-1, 1]))

            cur_predict_df = cur_predict_df.sort_values(by=[cur_celltype_column + '_TSG_prob'], ascending=False)
            cur_predict_df[cur_celltype_column + '_TSG_prob'] = np.arange(0, 1, 1. / cur_predict_df.shape[0])
            cur_predict_df[cur_celltype_column + '_TSG_prob'] = 1. / (1 + (
            0.05 / (cur_predict_df[cur_celltype_column + '_TSG_prob'] + eps))) ** 4
            cur_predict_df[cur_celltype_column + '_TSG_prob'] = 1 - MinMaxScaler().fit(
                cur_predict_df[cur_celltype_column + '_TSG_prob'].values.reshape([-1, 1])).transform(
                cur_predict_df[cur_celltype_column + '_TSG_prob'].values.reshape([-1, 1]))


            # df['rank'] = np.arange(0, 1, 1./df.shape[0])
            # df['rank'] = 1./(1+(0.05/(df['rank']+eps)))**4
            # df['rank'] = 1-MinMaxScaler().fit(df['rank'].values.reshape([-1,1])).transform(df['rank'].values.reshape([-1,1]))


            cur_predict_df['context_score'] = cur_context_df['distance']

            og_candidates = list(cur_predict_df.nlargest(500, cur_celltype_column+'_OG_prob').index)
            tsg_candidates = list(cur_predict_df.nlargest(500, cur_celltype_column+'_TSG_prob').index)

            # pan_og_candidates = list(cur_predict_df.nlargest(500, 'OG_prob').index)
            # pan_tsg_candidates = list(cur_predict_df.nlargest(500, 'TSG_prob').index)

            result_df.ix[cur_celltype_column, 'p_og'] = -np.log10(stats.mannwhitneyu(exp_df.ix[og_candidates, 'HMEL_BREAST'], exp_df.ix[og_candidates, cur_exp_column], alternative='less')[1])
            result_df.ix[cur_celltype_column, 'pan_p_og'] = -np.log10(stats.mannwhitneyu(exp_df.ix[pan_og_candidates, 'HMEL_BREAST'], exp_df.ix[pan_og_candidates, cur_exp_column],
                                    alternative='less')[1])

            result_df.ix[cur_celltype_column, 'p_tsg'] = -np.log10(stats.mannwhitneyu(exp_df.ix[tsg_candidates, 'HMEL_BREAST'], exp_df.ix[tsg_candidates, cur_exp_column], alternative='greater')[1])
            result_df.ix[cur_celltype_column, 'pan_p_tsg'] =  -np.log10(stats.mannwhitneyu(exp_df.ix[pan_tsg_candidates, 'HMEL_BREAST'], exp_df.ix[pan_tsg_candidates, cur_exp_column],
                                    alternative='greater')[1])


            # tsg_median = exp_df.ix[tsg_candidates, 'HMEL_BREAST'].median(), exp_df.ix[tsg_candidates, cur_exp_column].median()
            # og_median = exp_df.ix[og_candidates, 'HMEL_BREAST'].median(), exp_df.ix[og_candidates, cur_exp_column].median()
            #
            # pan_tsg_median = exp_df.ix[pan_tsg_candidates, 'HMEL_BREAST'].median(), exp_df.ix[pan_tsg_candidates, cur_exp_column].median()
            # pan_og_median = exp_df.ix[pan_og_candidates, 'HMEL_BREAST'].median(), exp_df.ix[pan_og_candidates, cur_exp_column].median()


            # for gene in tsg_candidates:
            #     print exp_df.ix[gene, :]
            cur_predict_df['HMEL_BREAST'] = exp_df['HMEL_BREAST']
            cur_predict_df[cur_exp_column] = exp_df[cur_exp_column]

            cur_predict_df.to_excel('../data/breast/'+cur_celltype_column+'_adjust_prediction.xlsx')

            cur_predict_df.nlargest(500, cur_celltype_column+'_OG_prob').to_excel('../data/breast/'+cur_celltype_column+'_adjust_OG_prediction.xlsx')

            cur_predict_df.nlargest(500, cur_celltype_column+'_TSG_prob').to_excel('../data/breast/'+cur_celltype_column+'_adjust_TSG_prediction.xlsx')

            og_df = known_gtf[known_gtf['hg19.kgXref.geneSymbol'].isin(og_candidates)]
            og_df.to_csv(cur_celltype_column+'_OGs.txt', sep='\t', header=None, index=None)

            tsg_df = known_gtf[known_gtf['hg19.kgXref.geneSymbol'].isin(tsg_candidates)]
            tsg_df.to_csv(cur_celltype_column + '_TSGs.txt', sep='\t', header=None, index=None)

            control_candidates = cur_predict_df[(~cur_predict_df.index.isin(og_candidates)) & (~cur_predict_df.index.isin(tsg_candidates)) & (~cur_predict_df.index.isin(pan_og_candidates)) & (~cur_predict_df.index.isin(pan_tsg_candidates))].sample(500, random_state=0).index

            control_df = known_gtf[known_gtf['hg19.kgXref.geneSymbol'].isin(control_candidates)]
            control_df.to_csv(cur_celltype_column + '_controls.txt', sep='\t', header=None, index=None)

        result_df.to_csv('../data/breast/p_values.csv')


    """
    common oncogenes and tumor suppressors based on different types of tumor celltypes
    """
    if True:
        celltypes =[['MCF7', 'ZR751'],
                    ['MB361', 'UACC812'],
                    ['SKBR3','AU565','HCC1954'],
                    ['MB468','HCC1937'],
                    ['MB231','MB436'],
                    ['MB468', 'HCC1937', 'MB231', 'MB436']]
        path = '../data/breast/'
        common_dict = {}
        for celltype in celltypes:
            name = '_'.join(celltype)
            common_dict[name] = {}
            common_dict[name]['OG'] = set()
            common_dict[name]['TSG'] = set()

            for c in celltype:
                cur_tsg_df = pd.read_excel(path+c+'_adjust_TSG_prediction.xlsx', index_col=0)
                cur_og_df = pd.read_excel(path + c + '_adjust_OG_prediction.xlsx', index_col=0)

                if len(common_dict[name]['TSG'])==0:
                    common_dict[name]['TSG'] = common_dict[name]['TSG'].union(cur_tsg_df.index)
                    common_dict[name]['OG'] = common_dict[name]['OG'].union(cur_og_df.index)
                else:
                    common_dict[name]['TSG'] = common_dict[name]['TSG'].intersection(cur_tsg_df.index)
                    common_dict[name]['OG'] = common_dict[name]['OG'].intersection(cur_og_df.index)
            common_dict[name]['TSG'] = list(common_dict[name]['TSG'])
            common_dict[name]['OG'] = list(common_dict[name]['OG'])

        for celltype in common_dict.keys():
            cur_shape = max(len(common_dict[celltype]['OG']), len(common_dict[celltype]['TSG']))
            cur_df = pd.DataFrame(index=range(cur_shape))
            cur_df['OG'] = common_dict[celltype]['OG'] + ['']*(cur_shape - len(common_dict[celltype]['OG']))
            cur_df['TSG'] = common_dict[celltype]['TSG'] + [''] * (cur_shape - len(common_dict[celltype]['TSG']))
            cur_df.to_excel(path+celltype+'_cancer_genes.xlsx', index=None)





