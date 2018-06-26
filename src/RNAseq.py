import pandas as pd, numpy as np
from scipy import stats

if __name__ == "__main__":
    """
    Extract the breast cancer cell line RPKM from CCLE to match the chipseq data from GSE96867
    """

    if True:
        df_path = '../data/CCLE_CancerCellLine_RNAseq_RPKM_20180214.xls'
        ccle_df = pd.read_csv(
            '/Users/boxia/PycharmProjects/CIG_logistic_regression/Fig8&Sup/CCLE_CancerCellLine_RNAseq_RPKM_20180214.xls',
            sep='\t', header=2, index_col=1)
        ccle_df = ccle_df[[x for x in ccle_df.columns if x.find('_BREAST')!=-1]]

        ccle_df = ccle_df[~ccle_df.index.duplicated(keep='first')]

        ccle_df.to_csv('../data/breast/CCLE_breast.xls', sep='\t')



