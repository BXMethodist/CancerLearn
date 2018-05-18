import pandas as pd
from sklearn import preprocessing

def center_normalization(table):
    return preprocessing.StandardScaler().fit(table)

def preprocessing_table(scaler, table):
    new_table_data = scaler.transform(table)
    new_table = pd.DataFrame(new_table_data, index=table.index, columns=table.columns)
    return new_table

def any_row(df):
    results = []
    for i in df.index:
        for c in df.columns:
            if df.ix[i, c] < 0.4:
                results.append(i)
                break
    return results


if __name__ == '__main__':
    ## Get the TSG and OG by using the PAN-cancer q-value and tissue specific p-values
    if False:
        tsg_df = pd.read_excel('../train/TSG.xlsx', index_col=0)
        # bq_tsg_df = tsg_df[(tsg_df['Breast_q-value']<0.05)]
        bp_tsg_df = tsg_df[(tsg_df['Breast_p-value']<0.05) | (tsg_df['PAN-Cancer_q-value']<0.05)]
        # bq_tsg_df[['Breast_q-value']].to_csv('TSG_breast_q.xls', sep='\t')
        bp_tsg_df[['Breast_p-value', 'PAN-Cancer_p-value']].to_csv('TSG_breast_p.xls', sep='\t')

        # pq_tsg_df = tsg_df[(tsg_df['Prostate_q-value'] < 0.05)]
        pp_tsg_df = tsg_df[(tsg_df['Prostate_p-value'] < 0.05) | ((tsg_df['PAN-Cancer_q-value']<0.05))]
        # pq_tsg_df[['Prostate_q-value']].to_csv('TSG_Prostate_q.xls', sep='\t')
        pp_tsg_df[['Prostate_p-value', 'PAN-Cancer_p-value']].to_csv('TSG_Prostate_p.xls', sep='\t')

        og_df = pd.read_excel('../train/OG.xlsx', index_col=0)
        # bq_og_df = og_df[(og_df['Breast_q-value']<0.05)]
        bp_og_df = og_df[(og_df['Breast_p-value']<0.05) | ((og_df['PAN-Cancer_q-value']<0.05))]
        # bq_og_df[['Breast_q-value']].to_csv('OG_breast_q.xls', sep='\t')
        bp_og_df[['Breast_p-value', 'PAN-Cancer_p-value']].to_csv('OG_breast_p.xls', sep='\t')

        # pq_og_df = og_df[(og_df['Prostate_q-value'] < 0.05)]
        pp_og_df = og_df[(og_df['Prostate_p-value'] < 0.05) | ((og_df['PAN-Cancer_q-value']<0.05))]
        # pq_og_df[['Prostate_q-value']].to_csv('OG_Prostate_q.xls', sep='\t')
        pp_og_df[['Prostate_p-value', 'PAN-Cancer_p-value']].to_csv('OG_Prostate_p.xls', sep='\t')


    ## Get the TSG and OG top 500 by using the minimum p values across different tissues
    if True:
        tsg_df = pd.read_excel('../genelist/TSG.xlsx', index_col=0)
        tsg_df = tsg_df[[x for x in tsg_df.columns if x.find('p-value')!=-1]]
        og_df = pd.read_excel('../genelist/OG.xlsx', index_col=0)
        og_df = og_df[[x for x in og_df.columns if x.find('p-value')!=-1]]

        control_df = pd.read_excel('../genelist/final_control_genes.xlsx', index_col=0)

        for i in range(500, 0, -100):
            print i
            tsg_candidates = tsg_df.nsmallest(i, 'PAN-Cancer_p-value').index

            og_df = og_df.sort_values(by='PAN-Cancer_p-value')
            og_candidates = og_df.ix[:i, :].index

            print len(tsg_candidates), len(og_candidates)
            overlap_tsg = any_row(tsg_df.ix[og_candidates, :])
            overlap_og = any_row(og_df.ix[tsg_candidates, :])

            og_df = og_df.ix[og_candidates, ['Breast_p-value', 'Prostate_p-value', 'PAN-Cancer_p-value', ]]
            tsg_df = tsg_df.ix[tsg_candidates, ['Breast_p-value', 'Prostate_p-value', 'PAN-Cancer_p-value', ]]

            og_df = og_df[~og_df.index.isin(overlap_tsg)]
            tsg_df = tsg_df[~tsg_df.index.isin(overlap_og)]

            og_df['type'] = ['OG']*og_df.shape[0]
            tsg_df['type'] = ['TSG']*tsg_df.shape[0]

            og_df['label'] = [2] * og_df.shape[0]
            tsg_df['label'] = [1] * tsg_df.shape[0]

            print og_df.shape
            print tsg_df.shape

            cancer_df = tsg_df.append(og_df).append(control_df)
            cancer_df.to_excel('../genelist/TUSON_top'+str(i)+'.xlsx')



    ## Get the control gene list by using the minimum p values across different tissues > 0.5
    if False:
        tsg_df = pd.read_excel('../genelist/TSG.xlsx', index_col=0)
        tsg_df = tsg_df[[x for x in tsg_df.columns if x.find('p-value')!=-1]]
        tsg_df['min'] = tsg_df.min(axis=1)
        candidates1 = tsg_df[tsg_df['min']>0.45].index
        # print len(candidates1)

        og_df = pd.read_excel('../genelist/OG.xlsx', index_col=0)
        og_df = og_df[[x for x in og_df.columns if x.find('p-value')!=-1]]
        og_df['min'] = og_df.min(axis=1)
        candidates2 = og_df[og_df['min']>0.4].index
        # print len(candidates2)

        candidates = set(candidates1).intersection(set(candidates2))
        # print len(candidates)
        results = []
        for candidate in candidates:
            results.append([candidate, 'Control', 0])
        df = pd.DataFrame(results)
        df.columns =['gene', 'type', 'label']
        df.to_excel('../genelist/Control_genes.xlsx', index=None)


    ## Get the TSG and OG from COSMIC project
    if False:
        df = pd.read_csv('../genelist/cancer_gene_census.csv', index_col=0)
        df = df.fillna('')

        og_df = df[(df['Role in Cancer'].str.contains('oncogene'))&(~(df['Role in Cancer'].str.contains('TSG'))) & (~(df['Role in Cancer'].str.contains('fusion')))]
        tsg_df = df[(df['Role in Cancer'].str.contains('TSG')) & (~(df['Role in Cancer'].str.contains('oncogene'))) & (~(df['Role in Cancer'].str.contains('fusion')))]

        og_df.to_csv('../genelist/cosmic_oncogene.xls', sep='\t')
        tsg_df.to_csv('../genelist/cosmic_tsg.xls', sep='\t')

        # og_prostate_df.to_csv('cosmic_prostate_oncogene.csv')
        # tsg_prostate_df.to_csv('cosmic_prostate_tsg.csv')



    ## Analyze the feature of the TSG and OG from COSMIC and TUSON
    if False:
        tuson_tsg = pd.read_csv('../genelist/TSG_top500.xls', sep='\t', index_col=0)
        tuson_og = pd.read_csv('../genelist/OG_top500.xls', sep='\t', index_col=0)

        default_df = pd.read_csv('../train/Cancer_CD4_default.csv', index_col=0)
        default_df.index = default_df.index.str.upper()
        cancer_df = default_df[default_df.index.isin(tuson_tsg.index)|default_df.index.isin(tuson_og.index)]

        columns = sorted([x for x in cancer_df.columns if x.find('genebody')!=-1 and x.find('width')!=-1 and x.find('single')==-1])
        cancer_df = cancer_df[columns]

        scaler = center_normalization(cancer_df)

        cancer_df = preprocessing_table(scaler, cancer_df)

        og_df = cancer_df[(cancer_df.index.isin(tuson_og.index)) & (~cancer_df.index.isin(tuson_tsg.index))]
        tsg_df = cancer_df[(cancer_df.index.isin(tuson_tsg.index))& (~cancer_df.index.isin(tuson_og.index))]

        print og_df.shape, tsg_df.shape

        overlap_df = cancer_df[(cancer_df.index.isin(tuson_tsg.index))& (cancer_df.index.isin(tuson_og.index))]

        cancer_df.T.to_excel('TUSON_Cancer_genes_default_genebody_width_values.xlsx')
        og_df.T.to_excel('TUSON_OG_genes_default_genebody_width_values.xlsx')
        tsg_df.T.to_excel('TUSON_TSG_genes_default_genebody_width_values.xlsx')
        overlap_df.T.to_excel('TUSON_overlap_default_genebody_with_values.xlsx')

        # control_df = pd.read_excel('../genelist/Control_genes.xlsx', index_col=0)
        # control_df.index = control_df.index.str.upper()
        # control_df = control_df.sample(100, random_state=0)