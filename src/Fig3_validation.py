import pandas as pd, numpy as np, os
from scipy import stats
import matplotlib.pyplot as plt
from utils import quantileNormalize

if __name__ == "__main__":
    """
    check whether the putative genes are enriched in TUSON databases
    Do the top overlap for bin size equal to 200
    """
    if False:
        predict_df = pd.read_excel('../results/aggregate_5_vote100.xlsx', index_col=0)
        TSG_df = predict_df.sort_values(by='TSG_prob', ascending=False).copy()
        OG_df = predict_df.sort_values(by='OG_prob', ascending=False).copy()

        tsg_df = pd.read_excel('../genelist/TSG.xlsx', index_col=0)
        tsg_df = tsg_df[[x for x in tsg_df.columns if x.find('p-value') != -1]]
        og_df = pd.read_excel('../genelist/OG.xlsx', index_col=0)
        og_df = og_df[[x for x in og_df.columns if x.find('p-value') != -1]]

        TUSON_TSG = tsg_df.nsmallest(200, 'PAN-Cancer_p-value').index

        og_df = og_df.sort_values(by='PAN-Cancer_p-value')
        TUSON_OG = og_df.ix[:200, :].index

        n = 0
        results = []
        for i in range(0, predict_df.shape[0], 200):
            start = i
            end = i + 200 if i+200 < predict_df.shape[0] else predict_df.shape[0]

            cur_TSG = TSG_df.iloc[start:end, :].index
            cur_OG = OG_df.iloc[start:end, :].index

            tsg_overlap = len(set(TUSON_TSG).intersection(cur_TSG))
            og_overlap = len(set(TUSON_OG).intersection(cur_OG))

            tsg_to_og_overlap = len(set(TUSON_OG).intersection(cur_TSG))
            og_to_tsg_overlap = len(set(TUSON_TSG).intersection(cur_OG))

            tsg_p = stats.fisher_exact([[tsg_overlap, 200-tsg_overlap],[200-tsg_overlap, predict_df.shape[0]-400+tsg_overlap]])[1]
            og_p = stats.fisher_exact(
                [[og_overlap, 200 - og_overlap], [200 - og_overlap, predict_df.shape[0] - 400 + og_overlap]])[1]
            tsg_to_og_p = stats.fisher_exact(
                [[tsg_to_og_overlap, 200 - tsg_to_og_overlap],
                 [len(TUSON_OG) - tsg_to_og_overlap, predict_df.shape[0] - 200 - len(TUSON_OG) + tsg_to_og_overlap]])[
                1]
            og_to_tsg_p = stats.fisher_exact(
                [[og_to_tsg_overlap, 200 - og_to_tsg_overlap], [len(TUSON_TSG) - og_to_tsg_overlap,
                                                                predict_df.shape[0] - 200 - len(
                                                                    TUSON_TSG) + og_to_tsg_overlap]])[1]
            results.append([n,
                            -np.log10(tsg_p),
                            -np.log10(og_p),
                            -np.log10(tsg_to_og_p),
                            -np.log10(og_to_tsg_p), ])
            n += 1

        final_df = pd.DataFrame(results)
        final_df.columns = ['bin', 'TSG_p', 'OG_p', 'TSG_to_OG', 'OG_to_TSG']
        final_df.to_excel('../plots/predict_TUSON_overlap.xlsx', index=None)

    """
        check whether the putative genes are enriched in COSMIC databases
        Do the top overlap for bin size equal to 200
    """
    if False:
        predict_df = pd.read_excel('../results/aggregate_5_vote100.xlsx', index_col=0)
        TSG_df = predict_df.sort_values(by='TSG_prob', ascending=False).copy()
        OG_df = predict_df.sort_values(by='OG_prob', ascending=False).copy()

        df = pd.read_csv('../genelist/cancer_gene_census.csv', index_col=0)
        df = df.fillna('')

        # og_df = df[
        #     (df['Role in Cancer'].str.contains('oncogene')) & (~(df['Role in Cancer'].str.contains('TSG'))) & (
        #     ~(df['Role in Cancer'].str.contains('fusion')))]
        # tsg_df = df[
        #     (df['Role in Cancer'].str.contains('TSG')) & (~(df['Role in Cancer'].str.contains('oncogene'))) & (
        #     ~(df['Role in Cancer'].str.contains('fusion')))]
        COSMIC_OG = df[
            (df['Role in Cancer'].str.contains('oncogene'))].index
        COSMIC_TSG = df[
            (df['Role in Cancer'].str.contains('TSG'))].index

        n = 0
        results = []
        for i in range(0, predict_df.shape[0], 200):
            start = i
            end = i + 200 if i + 200 < predict_df.shape[0] else predict_df.shape[0]

            cur_TSG = TSG_df.iloc[start:end, :].index
            cur_OG = OG_df.iloc[start:end, :].index

            tsg_overlap = len(set(COSMIC_TSG).intersection(cur_TSG))
            og_overlap = len(set(COSMIC_OG).intersection(cur_OG))

            tsg_to_og_overlap = len(set(COSMIC_OG).intersection(cur_TSG))
            og_to_tsg_overlap = len(set(COSMIC_TSG).intersection(cur_OG))

            tsg_p = stats.fisher_exact(
                [[tsg_overlap, 200 - tsg_overlap], [len(COSMIC_TSG) - tsg_overlap, predict_df.shape[0] - 200-len(COSMIC_TSG) + tsg_overlap]])[1]
            og_p = stats.fisher_exact(
                [[og_overlap, 200 - og_overlap], [len(COSMIC_OG) - og_overlap, predict_df.shape[0] - 200-len(COSMIC_OG) + og_overlap]])[1]

            tsg_to_og_p = stats.fisher_exact(
                [[tsg_to_og_overlap, 200 - tsg_to_og_overlap], [len(COSMIC_OG) - tsg_to_og_overlap, predict_df.shape[0] - 200-len(COSMIC_OG) + tsg_to_og_overlap]])[1]
            og_to_tsg_p = stats.fisher_exact(
                [[og_to_tsg_overlap, 200 - og_to_tsg_overlap], [len(COSMIC_TSG) - og_to_tsg_overlap, predict_df.shape[0] - 200-len(COSMIC_TSG) + og_to_tsg_overlap]])[1]
            results.append([n,
                            -np.log10(tsg_p),
                            -np.log10(og_p),
                            -np.log10(tsg_to_og_p),
                            -np.log10(og_to_tsg_p),])
            n += 1

        final_df = pd.DataFrame(results)
        final_df.columns = ['bin', 'TSG_p', 'OG_p', 'TSG_to_OG', 'OG_to_TSG']
        final_df.to_excel('../plots/predict_COSMIC_overlap.xlsx', index=None)







