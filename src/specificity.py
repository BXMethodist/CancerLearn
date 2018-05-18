import pandas as pd, numpy as np, os




tau_index = (1-final_df.iloc[:, 1:].divide(final_df.iloc[:, 1:].max(axis=1), axis=0)).sum(axis=1)/9