from .run_hit_func import *

def find_hits(ip, primary, verbose=True, correction=True):
    """ This function is used to take in pre-computed subcategorized data and calculate the chi-square contingency table 
    of a single gene or probe with all other genes or probes

    ip Input subcategorized data with 1 or 2s
    primary The probe or gene used to calculate chi-square contingency table with all other genes

    Returns
     - p-value of all matches
     - p-value <= 0.05 for all matches
     - table of the 2x2table, pvalue, rvalue

    """
    p_val = []
    odds_ratio = []
    ipa = ip.copy().T
    ipa.replace(3, 2, inplace=True)  # replace string with integer
    # after substituting some are then all 2s
    ipa = ipa.loc[:, (ipa != 2).any(axis=0)]

    for x in tqdm(ipa.columns):  # DONT DO ANYTHING WITH LOCALIZATION

        ipan = ipa[ipa[x] != np.nan]  # take out uncertain = 0

        try:
            o, p = stats.fisher_exact(pd.crosstab(ipan[primary], ipan[x]))
            p_val.append(p)
            odds_ratio.append(o)
            del o, p  # free memory
        except:
            p_val.append(1)
            odds_ratio.append(1)

    new = pd.DataFrame(
        {'P-value': p_val}, index=ipa.columns).sort_values('P-value', ascending=True)

    if correction == True:
        filtnew = new[(new < (0.05/len(ip)))['P-value']]
    else:
        filtnew = new[(new < 0.05)['P-value']]

    # only return those that are significant

    index = filtnew.index.drop(primary)

    ipa = ip.copy().T

    ct = []
    # take out uncertain = 0

    for x in index:
        ipan = ipa[~np.isnan(ipa[x])]
        o, p = stats.fisher_exact(pd.crosstab(ipan[primary], ipan[x]))

        first = ipan[primary].dropna()

        r, pp = stats.pearsonr(first, ipan.loc[first.index, x])

        # extract frm crosstab
        values = [y for x in pd.crosstab(
            ipan[primary], ipan[x]).values for y in x]
        ct.append([values[3], values[2], values[1], values[0],
                   p, r, inclusion_criterion(values, r)])
        if verbose == True:
            print(pd.crosstab(ipan[primary], ipan[x]))
            print('P-value: %s' % p+'\n')
            print('R-value: %s' % r+'\n')

    return new, filtnew, ct
