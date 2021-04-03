import numpy as np
import pandas as pd
import time
from tqdm import tqdm_notebook as tqdm
import math
from scipy import stats

# chisquare formula
def chi_square(o, e):
    return (o-e)**2/e

# assess whether it is an inclusion


def inclusion_criterion(ct, r):
    import numpy as np
    from scipy.stats import chi2
    equality = stats.chi2.ppf(0.95, 1)

    if r > 0:
        var1 = ct[1]
        var2 = ct[2]
        chi = (var1-var2)**2/sum([var1, var2])

        if chi > equality:
            return 'Significant'
        else:
            return 'No'
    else:
        var1 = ct[0]
        var2 = ct[3]
        chi = (var1-var2)**2/sum([var1, var2])

        if chi > equality:
            return 'Significant'
        else:
            return 'No'

# run all hits against one
