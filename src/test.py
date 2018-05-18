import pandas as pd

df = pd.read_csv('../ref_data/PolyA.1000.unstable.P.txt', header=None, index_col=0, sep='\t')
predict_df = pd.read_excel('../results/aggregated_real_results.xlsx', index_col=0)

predict_df=predict_df[predict_df.index.isin(df.index)]

print predict_df.to_string()