"""
Fig 1 to describe the schema of the model
"""

import pandas as pd, numpy as np, os

if __name__ == "__main__":
    """
    Get the variation presentation
    """
    if True:
        tables = [x for x in os.listdir('../feature_tables/') if x.endswith('.csv')]
        features =['H3K4me3_total_width', 'H3K27me3_total_width_genebody']

        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        OGs = cancer_df[(cancer_df['label'] == 2)].index
        TSGs = cancer_df[(cancer_df['label'] == 1)].index

        for feature in features:
            cur_table = '../feature_tables/'+feature+'.csv'
            cur_df = pd.read_csv(cur_table, index_col=0)

            cur_og_df = cur_df[cur_df.index.isin(OGs)]
            cur_tsg_df = cur_df[cur_df.index.isin(TSGs)]

            cur_og_df.to_csv('../plots/'+feature+'_OG.csv')
            cur_tsg_df.to_csv('../plots/' + feature + '_TSG.csv')


