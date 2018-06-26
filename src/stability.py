"""
Create the matrix for all genes, group by features and calculate CV, SD, tau index
"""
import pandas as pd, os, numpy as np
from utils import *
from sklearn import svm
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE

if __name__ == "__main__":
    """
    Generate tables for STD, CV, TAU, MIN, MAX, MEAN for each epigenetic feature from different cell types.
    """
    if False:
        path = '../real_tables/5/'
        real_tables = [x for x in os.listdir(path) if x.endswith('_real_table.csv')]

        # candidates_celltypes = ['HUVEC', 'HMEC', 'HSMM', 'NHLF']

        features = {}
        for table in real_tables:
            # if table.split('_')[0] not in candidates_celltypes:
            #     continue
            df = pd.read_csv(path+table, index_col=0)
            for c in df.columns:
                # if c.find('H3K9me3')!=-1:
                #     continue
                features[c] = pd.DataFrame(index=df.index)
            break
        for feature in features.keys():
            if feature.find('single')!=-1:
                continue
            cur_df = features[feature]
            for table in real_tables:
                print table
                # if table.split('_')[0] not in candidates_celltypes:
                #     continue
                celltype = table.replace('_parameters_default_real_table.csv', '')
                cur_table = pd.read_csv(path + table, index_col=0)
                cur_df[celltype+'_'+feature] = cur_table[feature]

            # cur_df = quantileNormalize(cur_df)
            tau_index = tau(cur_df)
            std = cur_df.std(axis=1)
            mean = cur_df.mean(axis=1)
            cv = std/mean
            max = cur_df.max(axis=1)
            min = cur_df.min(axis=1)
            cur_df['tau'] =tau_index
            cur_df['std'] = std
            cur_df['mean'] = mean
            cur_df['cv'] = cv
            cur_df['max'] = max
            cur_df['min'] = min
            cur_df.to_csv('../feature_tables/'+feature+'.csv')

    """
    add control genes to training data
    """
    if False:
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        cancer_df = cancer_df[(cancer_df['label']==1)|(cancer_df['label']==2)]


        all_genes = pd.read_csv('../real_tables/5/CD4_parameters_default_real_table.csv', index_col=0)
        all_genes = all_genes[~all_genes.index.isin(cancer_df)]
        controls = all_genes.sample(150, random_state=100)

        # controls = pd.read_excel('../genelist/final_control_genes.xlsx', index_col=0)
        controls_df = pd.DataFrame(index=controls.index)
        controls_df['type'] = ['Control']*controls_df.shape[0]
        controls_df['label'] = [0] * controls_df.shape[0]
        cancer_df = cancer_df.append(controls_df)
        cancer_df = cancer_df[~cancer_df.index.duplicated(keep='first')]
        cancer_df.to_excel('../genelist/TUSON_cancer_genes_curated.xlsx')

    """
    Check whether there is difference between tumor suppressor and oncogene and control
    """
    if False:
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)

        training_df = pd.DataFrame(index=cancer_df.index)
        og = cancer_df[cancer_df.type=='OG'].index
        tsg = cancer_df[cancer_df.type == 'TSG'].index
        control = cancer_df[cancer_df.type == 'Control'].index
        feature_tables = [x for x in os.listdir('../feature_tables') if x.endswith('.csv')]
        results = []
        for feature in feature_tables:
            print feature
            # if feature.find('gene_length')!=-1:
            #     continue
            cur_df = pd.read_csv('../feature_tables/'+feature, index_col=0)
            cur_og = cur_df[cur_df.index.isin(og)]
            cur_tsg = cur_df[cur_df.index.isin(tsg)]
            cur_control = cur_df[cur_df.index.isin(control)]
            cur_feature = feature[:-4]
            columns = ['tau','std','mean','cv','max','min',]
            for c in columns:
                p1 = logP_wilcoxon(cur_og[c], cur_tsg[c])
                p2 = logP_wilcoxon(cur_control[c], cur_tsg[c])
                p3 = logP_wilcoxon(cur_og[c], cur_control[c])
                if p1 < 0 or p2 <0 or p3<0:
                    results.append([cur_feature, c, p1, p2, p3])
                    training_df[feature + '_' + c] = cur_df[cur_df.index.isin(training_df.index)][c]
                    cur_result_df = pd.DataFrame(index=training_df.index)
                    cur_result_df['og'] = cur_og[c]
                    cur_result_df['tsg'] = cur_tsg[c]
                    cur_result_df['control'] = cur_control[c]
                    cur_result_df.to_excel('../results/'+cur_feature+'_'+c+'.xlsx')

                    ax = cur_result_df.boxplot(column=['tsg', 'og', 'control'])
                    fig = ax.get_figure()
                    fig.savefig('../results/'+cur_feature+'_'+c+'.pdf')
                    plt.close('all')

        result_df = pd.DataFrame(results)
        result_df.columns=['feature', 'agg', 'og_tsg', 'c_tsg', 'c_og']
        result_df.to_excel('../results/aggregated_features.xlsx', index=None)
        print training_df.columns
        training_df.to_excel('../results/aggregated_training.xlsx')

    """
       check the prediction result, without SMOTE
    """
    if False:
        from imblearn.over_sampling import SMOTE

        train_df = pd.read_excel('../results/aggregated_training.xlsx', index_col=0)
        train_df = train_df.fillna(0)
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top100.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top200.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top300.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top400.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top500.xlsx', index_col=0)
        train_df = label_label(train_df, cancer_df)

        scaler = center_normalization(train_df.iloc[:, :-1])

        train_df.iloc[:, :-1] = preprocessing_table(scaler, train_df.iloc[:, :-1])

        from CancerROC import ROC_plot
        from sklearn.linear_model import LogisticRegression

        # estimator = LogisticRegression(penalty='l1', C=.1)
        # estimator = svm.SVC(probability=True)
        estimator = RandomForestClassifier(n_estimators=10)
        model = OneVsRestClassifier(estimator)

        ROC_plot(model, train_df.iloc[:, :-1], train_df.iloc[:, -1], classes=[0, 1, 2], fname='wo_smote_ROC.pdf')

    """
    check the prediction result, with SMOTE, TOO GOOD TO BE TRUE
    """
    if False:
        train_df = pd.read_excel('../results/aggregated_training.xlsx', index_col=0)
        train_df = train_df.fillna(0)
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top100.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top200.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top300.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top400.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top500.xlsx', index_col=0)
        train_df = label_label(train_df, cancer_df)

        print train_df[train_df['label'] == 0].shape
        print train_df[train_df['label'] == 1].shape
        print train_df[train_df['label'] == 2].shape

        X = train_df.ix[:, :-1]
        Y = train_df.ix[:, -1]

        r2 = 6300
        X_resampled, Y_resampled = SMOTE(random_state=r2).fit_sample(X, Y)

        print X_resampled.shape

        smote_df = pd.DataFrame(X_resampled)
        smote_df.columns = train_df.columns[:-1]
        smote_df['label'] = Y_resampled
        smote_df.to_excel('../results/aggregated_training_smote.xlsx')

        print smote_df[smote_df['label']==0].shape
        print smote_df[smote_df['label'] == 1].shape
        print smote_df[smote_df['label'] == 2].shape

        scaler = center_normalization(smote_df.iloc[:, :-1])

        smote_df.iloc[:, :-1] = preprocessing_table(scaler, smote_df.iloc[:, :-1])

        from CancerROC import ROC_plot
        from sklearn.linear_model import LogisticRegression

        # estimator = LogisticRegression(penalty='l1', C=.1)
        # estimator = svm.SVC(probability=True)
        estimator = RandomForestClassifier(n_estimators=10)
        model = OneVsRestClassifier(estimator)

        ROC_plot(model, smote_df.iloc[:, :-1], smote_df.iloc[:, -1], classes=[0, 1, 2], fname='smote_ROC.pdf')

    """
    Check why smote can boost performance, no smote PCA analysis
    """
    if False:
        from mpl_toolkits.mplot3d import Axes3D
        from sklearn.decomposition import PCA

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        train_df = pd.read_excel('../results/aggregated_training.xlsx', index_col=0)
        train_df = train_df.fillna(0)
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)

        scaler = center_normalization(train_df)
        train_df = preprocessing_table(scaler, train_df)

        pca = PCA()
        pca.fit(train_df)
        print pca.explained_variance_ratio_
        train_df_new = train_df.copy()
        train_df_new.ix[:, :] = pca.fit_transform(train_df)
        # print train_df_new.index
        # pass

        train_df = label_label(train_df, cancer_df)
        # train_df1 = train_df1.sample(200, random_state=0)

        tsg_df = train_df_new[train_df.label==1]
        control_df = train_df_new[train_df.label==0]
        og_df = train_df_new[train_df.label == 2]

        # ax.scatter(control_df.ix[:, 0],
        #            control_df.ix[:, 1],
        #            control_df.ix[:, 2],
        #            c='grey',
        #            marker='s')
        ax.scatter(tsg_df.ix[:, 0],
                   tsg_df.ix[:, 1],
                   tsg_df.ix[:, 2],
                   c='red',
                   marker='o')
        ax.scatter(og_df.ix[:, 0],
                   og_df.ix[:, 1],
                   og_df.ix[:, 2],
                   c='blue',
                   marker='^')

        ax.set_xlabel('PCA1')
        ax.set_ylabel('PCA2')
        ax.set_zlabel('PCA3')


        plt.show()


    """
    Check why smote can boost performance, smote PCA analysis
    """
    if False:
        from mpl_toolkits.mplot3d import Axes3D
        from sklearn.decomposition import PCA

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        train_df = pd.read_excel('../results/aggregated_training_smote.xlsx', index_col=0)
        train_df = train_df.fillna(0)

        scaler = center_normalization(train_df.iloc[:, :-1])
        train_df.iloc[:, :-1] = preprocessing_table(scaler, train_df.iloc[:, :-1])

        pca = PCA()
        pca.fit(train_df.iloc[:, :-1])
        print pca.explained_variance_ratio_
        train_df_new = train_df.copy()
        train_df_new.ix[:, :-1] = pca.fit_transform(train_df.iloc[:, :-1])

        tsg_df = train_df_new[train_df.label==1]
        control_df = train_df_new[train_df.label==0]
        og_df = train_df_new[train_df.label == 2]

        # ax.scatter(control_df.ix[:, 0],
        #            control_df.ix[:, 1],
        #            control_df.ix[:, 2],
        #            c='grey',
        #            marker='s')
        ax.scatter(tsg_df.ix[:, 0],
                   tsg_df.ix[:, 1],
                   tsg_df.ix[:, 2],
                   c='red',
                   marker='o')
        ax.scatter(og_df.ix[:, 0],
                   og_df.ix[:, 1],
                   og_df.ix[:, 2],
                   c='blue',
                   marker='^')

        ax.set_xlabel('PCA1')
        ax.set_ylabel('PCA2')
        ax.set_zlabel('PCA3')


        plt.show()

    """
    Get all genes real table
    """
    if False:
        feature_tables = [x for x in os.listdir('../feature_tables') if x.endswith('.csv')]
        training_df = pd.DataFrame(index=pd.read_csv('../feature_tables/' + feature_tables[0], index_col=0).index)

        for feature in feature_tables:
            print feature
            cur_df = pd.read_csv('../feature_tables/' + feature, index_col=0)
            cur_feature = feature[:-4]
            columns = ['tau', 'std', 'mean', 'cv', 'max', 'min', ]
            for c in columns:
                training_df[feature + '_' + c] = cur_df[c]
        training_df.to_excel('../results/aggregated_real_5.xlsx')


    """
        add control genes to training data, create different smote training set based on different random seed
    """
    if False:
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        cancer_df = cancer_df[(cancer_df['label'] == 1) | (cancer_df['label'] == 2)]

        real_table = pd.read_excel('../results/aggregated_real_5.xlsx', index_col=0)
        real_table = real_table.fillna(0)

        for i in range(0, 10000, 100):
            all_genes = pd.read_csv('../real_tables/5/CD4_parameters_default_real_table.csv', index_col=0)
            all_genes = all_genes[~all_genes.index.isin(cancer_df)]
            controls = all_genes.sample(150, random_state=i)

            # controls = pd.read_excel('../genelist/final_control_genes.xlsx', index_col=0)
            controls_df = pd.DataFrame(index=controls.index)
            controls_df['type'] = ['Control'] * controls_df.shape[0]
            controls_df['label'] = [0] * controls_df.shape[0]
            cur_cancer_df = cancer_df.append(controls_df)
            cur_cancer_df = cur_cancer_df[~cur_cancer_df.index.duplicated(keep='first')]

            cur_table = real_table[real_table.index.isin(cur_cancer_df.index)]
            cur_table['label'] = cur_cancer_df['label']

            X = cur_table.ix[:, :-1]
            Y = cur_table.ix[:, -1]

            r2 = 6300
            X_resampled, Y_resampled = SMOTE(random_state=r2).fit_sample(X, Y)

            print X_resampled.shape

            smote_df = pd.DataFrame(X_resampled)
            smote_df.columns = cur_table.columns[:-1]
            smote_df['label'] = Y_resampled
            smote_df.to_excel('../results/aggregated_training_5_smote'+str(i)+'.xlsx')

    """
        Do the prediction
    """
    if False:
        real_table = pd.read_excel('../results/aggregated_real_5.xlsx', index_col=0)
        real_table = real_table.fillna(0)
        great_table = pd.read_csv('../ref_data/hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t',
                                  dtype={"hg19.kgXref.geneSymbol": str})
        great_table = great_table.set_index(['hg19.kgXref.geneSymbol'])
        real_table = real_table[real_table.index.isin(great_table.index)]

        # y_score = 0
        og_y_prob = pd.DataFrame(index=real_table.index)
        tsg_y_prob = pd.DataFrame(index=real_table.index)

        times = range(0, 10000, 100)

        for i in times:
            smote_df = pd.read_excel('../results/aggregated_training_5_smote'+str(i)+'.xlsx', index_col=0)

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
            estimator = RandomForestClassifier(n_estimators=10)
            model = OneVsRestClassifier(estimator)

            predictor = model.fit(smote_df.iloc[:, :-1], smote_df.iloc[:, -1])

            cur_real_table = preprocessing_table(scaler, real_table)

            # y_score += predictor.decision_function(cur_real_table)
            y_score = predictor.predict_proba(cur_real_table)
            og_y_prob['predict'+str(i)] = y_score[:, 2]
            tsg_y_prob['predict' + str(i)] = y_score[:, 1]

        og_y_prob.to_excel('OG_aggregate_5_vote100.xlsx')
        tsg_y_prob.to_excel('TSG_aggregate_5_vote100.xlsx')

    """
    check the cutoff of the probablility
    """
    if False:
        og_df = pd.read_excel('../results/OG_aggregate_5_vote100.xlsx', index_col=0)
        tsg_df = pd.read_excel('../results/TSG_aggregate_5_vote100.xlsx', index_col=0)

        control_df = 1 - og_df - tsg_df

        result_df=pd.DataFrame(index=og_df.index)
        result_df['TSG_prob'] = tsg_df.mean(axis=1)
        result_df['OG_prob'] = og_df.mean(axis=1)
        result_df['Control_prob'] = control_df.mean(axis=1)

        # for index in og_df.index:
        #     for c in og_df.columns:
        #         og_p, tsg_p, c_p = og_df.ix[index, c], tsg_df.ix[index, c], control_df.ix[index, c]
        #         if og_p > tsg_p and og_p > c_p:
        #             result_df.ix[index, 'OG_prob'] += 1
        #         elif og_p < tsg_p and tsg_p > c_p:
        #             result_df.ix[index, 'TSG_prob'] += 1
        #
        # result_df = result_df/100.
        result_df.to_excel('aggregate_5_vote100.xlsx')


    """
        check the prediction result, with SMOTE, For each biological markers
    """
    if True:
        markers = ['H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac', 'H3K79me2', 'H3K9me3',
                   'CTCF']

        train_df = pd.read_excel('../results/aggregated_training.xlsx', index_col=0)
        train_df = train_df.fillna(0)
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top100.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top200.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top300.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top400.xlsx', index_col=0)
        # cancer_df = pd.read_excel('../genelist/TUSON_top500.xlsx', index_col=0)
        train_df = label_label(train_df, cancer_df)

        scaler = center_normalization(train_df.iloc[:, :-1])

        train_df.iloc[:, :-1] = preprocessing_table(scaler, train_df.iloc[:, :-1])

        from CancerROC import ROC_plot

        for marker in markers:
            cur_train_df = train_df[[x for x in train_df.columns if x.find(marker)!=-1]]

            estimator = RandomForestClassifier(n_estimators=10)
            model = OneVsRestClassifier(estimator)

            ROC_plot(model, cur_train_df, train_df.iloc[:, -1], classes=[0, 1, 2], fname=marker+'_ROC.pdf')




