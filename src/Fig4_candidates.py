import pandas as pd, numpy as np, os
from scipy import stats

if __name__ == "__main__":
    """
    pipeline to select candidates
    1. remove genes from TUSON list
    2. get their expression in breast cancer cell lines
    3. get their survivals
    4. get their mutations
    5. get their mutations survival
    6. check reference
    """
    # Step 1
    if False:
        cancer_df = pd.read_excel('../genelist/TUSON_cancer_genes_curated.xlsx', index_col=0)
        print cancer_df.shape
        predict_df = pd.read_excel('../results/aggregate_5_vote100.xlsx', index_col=0)
        predict_df = predict_df[~predict_df.index.isin(cancer_df.index)]
        OGs = predict_df[predict_df['OG_prob']>0.7]
        TSGs = predict_df[predict_df['TSG_prob']>0.7]

        tsg_df = pd.read_excel('../genelist/TSG.xlsx', index_col=0)
        og_df = pd.read_excel('../genelist/OG.xlsx', index_col=0)

        TUSON_TSG = tsg_df[tsg_df['PAN-Cancer_q-value'] < 0.5].index

        TUSON_OG = og_df[og_df['PAN-Cancer_q-value'] < 0.5].index

        print OGs.shape, TSGs.shape

        OGs = OGs[(~OGs.index.isin(TUSON_OG)) & (~OGs.index.isin(TUSON_TSG))]
        TSGs = TSGs[(~TSGs.index.isin(TUSON_OG)) & (~TSGs.index.isin(TUSON_TSG))]

        print OGs.shape, TSGs.shape

        OGs.to_excel('../results/aggregate_5_vote100_OGs_notin_TUSON.xlsx')
        TSGs.to_excel('../results/aggregate_5_vote100_TSGs_notin_TUSON.xlsx')

    # Step 2:
    if False:
        exp_df = pd.read_excel('../celllines_rnaseq/GSE96860_FPKM.xlsx', index_col=0)

        OGs = pd.read_excel('../results/aggregate_5_vote100_OGs_notin_TUSON.xlsx', index_col=0)
        TSGs = pd.read_excel('../results/aggregate_5_vote100_TSGs_notin_TUSON.xlsx', index_col=0)

        OG_exp_df = exp_df[exp_df.index.isin(OGs.index)]
        TSG_exp_df = exp_df[exp_df.index.isin(TSGs.index)]
        print exp_df.columns
        median = exp_df.median(axis=0)
        MCF10A = exp_df['MCF10A_FPKM'].median()
        MDB231 = exp_df['MDA-MB-231_FPKM'].median()

        print OG_exp_df.shape, TSG_exp_df.shape

        OG_exp_df = OG_exp_df[OG_exp_df['MDA-MB-231_FPKM']>MDB231]
        TSG_exp_df = TSG_exp_df[TSG_exp_df['MCF10A_FPKM']>MCF10A]

        print OG_exp_df.shape, TSG_exp_df.shape

        OG_exp_df.to_excel('../results/aggregate_5_vote100_OGs_notin_TUSON_breast.xlsx')
        TSG_exp_df.to_excel('../results/aggregate_5_vote100_TSGs_notin_TUSON_breast.xlsx')

    # Check mutation and CNA for breast candidates
    if False:
        path = '/Users/boxia/PycharmProjects/tools/results/'
        breast_CNA = pd.read_csv(path+'breast_all_CNA.xls', sep='\t', index_col=0)
        breast_mutation = pd.read_csv(path+'breast_all_mutation.xls', sep='\t', index_col=0)

        breast_og_candidates = pd.read_excel('../results/aggregate_5_vote100_OGs_notin_TUSON_breast.xlsx', index_col=0).index
        breast_tsg_candidates = pd.read_excel('../results/aggregate_5_vote100_TSGs_notin_TUSON_breast.xlsx', index_col=0).index

        og_CNA = breast_CNA[breast_CNA.index.isin(breast_og_candidates)]
        tsg_CNA = breast_CNA[breast_CNA.index.isin(breast_tsg_candidates)]

        og_mutation = breast_mutation[breast_mutation.index.isin(breast_og_candidates)]
        tsg_mutation = breast_mutation[breast_mutation.index.isin(breast_tsg_candidates)]

        og_CNA.to_excel('../results/breast_OG_CNA.xlsx')
        tsg_CNA.to_excel('../results/breast_TSG_CNA.xlsx')

        og_mutation.to_excel('../results/breast_OG_mutation.xlsx')
        tsg_mutation.to_excel('../results/breast_TSG_mutation.xlsx')

    ## compare the mutation and CNA rates with control genes
    if False:
        tsg_df = pd.read_excel('../genelist/TSG.xlsx', index_col=0)
        og_df = pd.read_excel('../genelist/OG.xlsx', index_col=0)
        TUSON_TSG = tsg_df[tsg_df['PAN-Cancer_q-value'] < 0.5].index
        TUSON_OG = og_df[og_df['PAN-Cancer_q-value'] < 0.5].index


        path = '/Users/boxia/PycharmProjects/tools/results/'
        breast_CNA = pd.read_csv(path + 'breast_all_CNA.xls', sep='\t', index_col=0)
        breast_mutation = pd.read_csv(path + 'breast_all_mutation.xls', sep='\t', index_col=0)
        breast_mutation['SUM'] = breast_mutation.sum(axis=1)

        breast_CNA['amp'] = breast_CNA['amp_high'] + breast_CNA['amp_low']
        breast_CNA['del'] = breast_CNA['del_high']  + breast_CNA['del_low']

        breast_og_candidates = pd.read_excel('../results/aggregate_5_vote100_OGs_notin_TUSON.xlsx',
                                             index_col=0).index
        breast_tsg_candidates = pd.read_excel('../results/aggregate_5_vote100_TSGs_notin_TUSON.xlsx',
                                              index_col=0).index

        og_CNA = breast_CNA[breast_CNA.index.isin(breast_og_candidates)]
        tsg_CNA = breast_CNA[breast_CNA.index.isin(breast_tsg_candidates)]
        # control_CNA = breast_CNA[
        #     (~breast_CNA.index.isin(breast_og_candidates)) & (~breast_CNA.index.isin(breast_tsg_candidates))]
        control_CNA = breast_CNA[(~breast_CNA.index.isin(breast_og_candidates)) & (~breast_CNA.index.isin(breast_tsg_candidates))&(~breast_CNA.index.isin(TUSON_TSG))&(~breast_CNA.index.isin(TUSON_OG))]

        og_mutation = breast_mutation[breast_mutation.index.isin(breast_og_candidates)]
        tsg_mutation = breast_mutation[breast_mutation.index.isin(breast_tsg_candidates)]
        # control_mutation = breast_mutation[(~breast_mutation.index.isin(breast_og_candidates)) & (
        # ~breast_mutation.index.isin(breast_tsg_candidates))]
        # print control_mutation.shape
        control_mutation = breast_mutation[(~breast_mutation.index.isin(breast_og_candidates)) & (~breast_mutation.index.isin(breast_tsg_candidates))&(~breast_mutation.index.isin(TUSON_TSG))&(~breast_mutation.index.isin(TUSON_OG))]
        # print control_mutation.shape, len(TUSON_TSG), len(TUSON_OG)

        # cna_rows = control_CNA.shape[0]
        cna_rows = max(og_CNA.shape[0], tsg_CNA.shape[0])
        control_CNA = control_CNA.sample(cna_rows, random_state=10)
        
        amp_cna_df = pd.DataFrame(index=range(cna_rows))
        og_cna = list(og_CNA['amp']) + [np.nan]*(cna_rows-og_CNA.shape[0])
        tsg_cna = list(tsg_CNA['amp']) + [np.nan] * (cna_rows - tsg_CNA.shape[0])
        control_cna = list(control_CNA['amp'])
        amp_cna_df['OG'] = og_cna
        amp_cna_df['TSG'] = tsg_cna
        amp_cna_df['control'] = control_cna
        amp_cna_df.to_csv('../results/breast_og_tsg_control_amp.xls', sep='\t')

        # print amp_cna_df.mean(axis=0)
        print stats.mannwhitneyu(amp_cna_df['OG'], amp_cna_df['TSG'], alternative='greater')[1]
        print stats.mannwhitneyu(amp_cna_df['OG'], amp_cna_df['control'], alternative='greater')[1]


        del_cna_df = pd.DataFrame(index=range(cna_rows))
        og_cna = list(og_CNA['del']) + [np.nan] * (cna_rows - og_CNA.shape[0])
        tsg_cna = list(tsg_CNA['del']) + [np.nan] * (cna_rows - tsg_CNA.shape[0])
        control_cna = list(control_CNA['del'])
        del_cna_df['OG'] = og_cna
        del_cna_df['TSG'] = tsg_cna
        del_cna_df['control'] = control_cna
        del_cna_df.to_csv('../results/breast_og_tsg_control_del.xls', sep='\t')

        # print del_cna_df.mean(axis=0)
        print stats.mannwhitneyu(del_cna_df['TSG'], del_cna_df['OG'], alternative='greater')[1]
        print stats.mannwhitneyu(del_cna_df['TSG'], del_cna_df['control'], alternative='greater')[1]

        # mutation_rows = control_mutation.shape[0]
        mutation_rows = max(og_mutation.shape[0], tsg_mutation.shape[0])
        control_mutation = control_mutation.sample(mutation_rows, random_state=10)
        SUM_mutation_df = pd.DataFrame(index=range(mutation_rows))
        og_mutation = list(og_mutation['SUM']) + [np.nan] * (mutation_rows - og_mutation.shape[0])
        tsg_mutation = list(tsg_mutation['SUM']) + [np.nan] * (mutation_rows - tsg_mutation.shape[0])
        control_mutation = list(control_mutation['SUM'])
        SUM_mutation_df['OG'] = og_mutation
        SUM_mutation_df['TSG'] = tsg_mutation
        SUM_mutation_df['control'] = control_mutation
        SUM_mutation_df.to_csv('../results/breast_og_tsg_control_SUM.xls', sep='\t')

        # print SUM_mutation_df.mean(axis=0)
        print stats.mannwhitneyu(SUM_mutation_df['TSG'], SUM_mutation_df['control'], alternative='greater')[1]
        print stats.mannwhitneyu(SUM_mutation_df['OG'], SUM_mutation_df['control'], alternative='greater')[1]

    ## Get the candidate list for breast cancer, expression in cell line, survival, prediction score > 0.7
    if False:
        breast_og_candidates = pd.read_excel('../results/aggregate_5_vote100_OGs_notin_TUSON_breast.xlsx',
                                             index_col=0).index
        breast_tsg_candidates = pd.read_excel('../results/aggregate_5_vote100_TSGs_notin_TUSON_breast.xlsx',
                                              index_col=0).index

        breast_survival_df = pd.read_csv('/Users/boxia/PycharmProjects/tools/survival/breast/breast_survival_candidates.xls',
                                  sep='\t', index_col=0)
        og_breast_survival = breast_survival_df[breast_survival_df['type']=='OG candidate'].index
        tsg_breast_survival = breast_survival_df[breast_survival_df['type'] == 'TSG candidate'].index

        og_candidates = set(breast_og_candidates).intersection(og_breast_survival)
        tsg_candidates = set(breast_tsg_candidates).intersection(tsg_breast_survival)

        ## Get breast candidates expression in breast cell lines
        exp_df = pd.read_excel('../celllines_rnaseq/GSE96860_FPKM.xlsx', index_col=0)
        og_candidates = list(og_candidates)
        tsg_candidates = list(tsg_candidates)

        exp_df= exp_df.ix[og_candidates+tsg_candidates, ['MCF10A_FPKM', 'MDA-MB-231_FPKM']]

        for gene in og_candidates:
            if exp_df.ix[gene, 'MCF10A_FPKM'] < exp_df.ix[gene, 'MDA-MB-231_FPKM']:
                print gene, 'OG'

        for gene in tsg_candidates:
            if exp_df.ix[gene, 'MCF10A_FPKM'] > exp_df.ix[gene, 'MDA-MB-231_FPKM']:
                print gene, 'TSG'

        exp_df.to_csv('breast_candidates_FPKM.xls', sep='\t')






        












    # Check survival and best combination


    # Check the reference
