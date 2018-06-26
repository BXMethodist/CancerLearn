"""
Fig 2 explanation of the model
"""

import pandas as pd, numpy as np, os
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
from utils import *
from collections import defaultdict

if __name__ == "__main__":
    """
    Get the features contributions for the model
    """
    if False:
        smote_df = pd.read_excel('../results/aggregated_training_5_smote' + str(0) + '.xlsx', index_col=0)

        ###
        # smote_df = pd.read_excel('../results/aggregated_training.xlsx', index_col=0)
        # smote_df = smote_df.fillna(0)
        # cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        # smote_df = label_label(smote_df, cancer_df)
        ###


        scaler = center_normalization(smote_df.iloc[:, :-1])
        smote_df.iloc[:, :-1] = preprocessing_table(scaler, smote_df.iloc[:, :-1])

        from CancerROC import ROC_plot
        from sklearn.linear_model import LogisticRegression

        # estimator = LogisticRegression(penalty='l1', C=.1)
        # estimator = svm.SVC(probability=True)
        columns = [x.replace('.csv', '') for x in smote_df.columns[:-1]]
        result_df = pd.DataFrame(index=columns)

        for r in range(0, 10000, 100):
            estimator = RandomForestClassifier(n_estimators=10, random_state=r)
            model = OneVsRestClassifier(estimator)

            predictor = model.fit(smote_df.iloc[:, :-1], smote_df.iloc[:, -1])

            if 'predictor_1' not in result_df.columns:
                result_df['predictor_1'] = predictor.estimators_[0].feature_importances_
                result_df['predictor_2'] = predictor.estimators_[1].feature_importances_
                result_df['predictor_3'] = predictor.estimators_[2].feature_importances_
            else:
                result_df['predictor_1'] += predictor.estimators_[0].feature_importances_
                result_df['predictor_2'] += predictor.estimators_[1].feature_importances_
                result_df['predictor_3'] += predictor.estimators_[2].feature_importances_

        result_df = result_df/100.
        result_df['sum'] = result_df.sum(axis=1)
        result_df = result_df.sort_values(by='sum', ascending=False)
        result_df.to_excel('../plots/features_contributions.xlsx')

        contributions = {'predictor_1': defaultdict(float),
                         'predictor_2': defaultdict(float),
                         'predictor_3': defaultdict(float),
                         'sum': defaultdict(float),}
        for index in result_df.index:
            marker = index.split('_')[0]
            contributions['sum'][marker] += result_df.ix[index, 'sum']
            contributions['predictor_1'][marker] += result_df.ix[index, 'predictor_1']
            contributions['predictor_2'][marker] += result_df.ix[index, 'predictor_2']
            contributions['predictor_3'][marker] += result_df.ix[index, 'predictor_3']

        contributions_df = pd.DataFrame(contributions)
        contributions_df.to_excel('../plots/features_contributions.xlsx')


