import numpy as np
import pandas as pd
from tqdm import tqdm_notebook as tqdm
import math


def probe_filter(datainputnorm, log2transform=False, filt=0, threshold_filter=0.01, variance_filter=0.0125):

    if log2transform == True:
        # pseudocount so that we do not ignore gene with 0 expression
        ip_pseudo = datainputnorm.replace(0, np.nan)
        datainput = np.log2(ip_pseudo)

        data = datainput.copy()

        aa = []
        var = []
        for a in tqdm(range(0, len(data))):
            total = sum(data.iloc[a, :] >= filt)
            aa.append(total)

        aap = pd.DataFrame({'threshold': aa}, index=data.index)
        # remove samples with all sample below expression threshold
        aap = aap[aap['threshold'] > math.ceil(
            threshold_filter*len(datainput.T))]

        datainput = datainput.loc[aap.index]
        datainputnorm = datainputnorm.loc[aap.index]

        for a in tqdm(range(0, len(datainput))):
            varss = np.var(datainput.iloc[a, :])
            var.append(varss)

        varp = pd.DataFrame(
            {'variance': var, 'threshold': aap['threshold']}, index=datainput.index)
        varp = varp[(varp['variance'] < (max(var)-min(var))*variance_filter) & (varp['threshold']
                                                                                <= math.ceil(threshold_filter*len(datainput.T))+1)]  # remove sample with low variance

        finalpd = pd.concat([varp, aap], axis=1, join='inner')

        # remove sample with low variance and below threshold filter
        data = datainputnorm.drop(finalpd.index)
    else:
        data = datainputnorm.copy()
        data.replace(np.nan, filt-1, inplace=True)
        aa = []
        var = []
        for a in tqdm(range(0, len(data))):
            percent = sum(data.iloc[a, :] >= filt)
            aa.append(percent)

        aap = pd.DataFrame({'threshold': aa}, index=data.index)
        # remove samples with all sample below expression threshold
        aap = aap[aap['threshold'] > math.ceil(
            threshold_filter*len(datainputnorm.T))]

        datainputnorm = datainputnorm.loc[aap.index]

        for a in tqdm(range(0, len(datainputnorm))):
            varss = np.var(datainputnorm.iloc[a, :])
            var.append(varss)

        varp = pd.DataFrame(
            {'variance': var, 'threshold': aap['threshold']}, index=datainputnorm.index)
        varp = varp[(varp['variance'] < (max(var)-min(var))*variance_filter) & (varp['threshold']
                                                                                <= math.ceil(threshold_filter*len(datainputnorm.T))+1)]  # remove sample with low variance

        finalpd = pd.concat([varp, aap], axis=1, join='inner')

        # remove sample with low variance and below threshold filter
        data = datainputnorm.drop(finalpd.index)
    return data
