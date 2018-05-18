import pandas as pd
from collections import defaultdict

df = pd.read_excel('../grid/best_parameters.xlsx')

results = defaultdict(float)
results_para = {}
final = []

for i in range(df.shape[0]):
    marker, feature, genebody, target, upstream, downstream, height, logP = df.ix[i, :]
    if logP < results[(marker, feature, genebody)]:
        results[(marker, feature, genebody)] = logP
        results_para[(marker, feature, genebody)] = [upstream, downstream, height, logP]

for key, value in results_para.items():
    if key[2] =='TSS':
        final.append([key[0], key[1]] + value)
    else:
        final.append([key[0], '_'.join([key[1], key[2]])]+value)

para_df = pd.DataFrame(final)
para_df.columns = ['marker','feature','upstream','downstream','height','logP']
para_df.to_csv('../ref_data/best_parameters_grid.csv', index=None)


