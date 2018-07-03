import pandas as pd
import numpy as np
from scipy.stats import logistic
eps = np.finfo(float).eps
m=0.5
for i in np.arange(0, 2, 0.1):
    print i, '\t', 1./(1+(m/(i+eps)))**4