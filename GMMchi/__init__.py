# GMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM###########################
import matplotlib
from numpy import *
import itertools as itert
from scipy.stats import norm
from sklearn.mixture import GaussianMixture as GMM
import time
from tqdm import tqdm_notebook as tqdm
import math
import matplotlib.pyplot as plt
import sklearn.metrics as sklm
import sklearn.feature_selection as sklf
import sklearn.preprocessing as sklp
import warnings
from scipy.stats import chi2
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.ticker as tick
import matplotlib.patches as mpatches
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib as plt
import seaborn as sns
sns.set(color_codes=True)

class GMMchi():
    def solve_gaussians(m1, m2, cov1, cov2, s1, s2, minimum, maximum):
        std1 = np.sqrt(cov1)
        std2 = np.sqrt(cov2)

        a = 1/(2*std1**2) - 1/(2*std2**2)
        b = m2/(std2**2) - m1/(std1**2)
        c = m1**2 / (2*std1**2) - m2**2 / (2*std2**2) - np.log((std2*s1)/(std1*s2))
        # GET THE INTERSECTION THAT IS BETWEEN THE TWO MEANS
        return np.sort([x for x in np.roots([a[0], b[0], c[0]]) if x < maximum and x > minimum])


    def dynamic_binning(observed, binedges, threshold=5, filt=0, final=True):
        continued = 1
        x = 0  # dynamically update range
        orgi = observed.copy()

        while x < len(observed):
            try:
                restart = True  # loop until larger than threshold
                while restart:

                    if observed[x] < threshold:
                        observed[x] = observed[x]+observed[x+1]
                        observed = np.delete(observed, x+1)
                        binedges = np.delete(binedges, x+1)
                    else:
                        restart = False
                x += 1
            except:  # sometimes you will get bins with >threshold in the very last bin
                if observed[-1] < threshold:

                    observed[-1] = observed[-1]+observed[-2]
                    observed = np.delete(observed, -2)
                    binedges = np.delete(binedges, -2)

        if len(orgi) == len(observed):
            continued = 0

        if final == False:
            largedata = np.arange(len(observed))[np.in1d(
                observed, orgi)]  # get overlapping index

            try:
                # check whether it is spanning the entire dataset
                if min(largedata) != 0 or max(largedata)+1 != len(observed):

                    # expand the space where there is no tail problem to the number of bins using mann and wald and leave the rest
                    if len(largedata) == 1:
                        # mann and wald
                        nums = int(1.88*observed[largedata[0]]**(2/5))
                        x = largedata[0]
                        newbins = np.linspace(
                            binedges[largedata[0]], binedges[largedata[0]+1], num=nums)
                        binedges = np.delete(
                            binedges, (largedata[0], largedata[0]+2))
                        binedges = np.sort(np.concatenate((binedges, newbins)))
                    else:
                        # mann and wald
                        nums = int(
                            1.88*sum(observed[min(largedata):max(largedata)+1])**(2/5))
                        newbins = np.linspace(
                            binedges[min(largedata)], binedges[max(largedata)+1], num=nums)
                        binedges = np.delete(binedges, np.arange(
                            min(largedata), max(largedata)+2))
                        binedges = np.sort(np.concatenate((binedges, newbins)))
            except:
                pass  # if they are the same during the middle stage

        return observed, binedges, continued


    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    # v
    ##############################################################################################################################
    ##############################################################################################################################


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

    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    # v
    ##############################################################################################################################
    ##############################################################################################################################


    mpl.rcParams['figure.dpi'] = 300
    warnings.filterwarnings("ignore")

    %matplotlib inline


    def GMM_pipeline(data, log2transform,  filt, meanf, stdf, verbosity=False, farinto=0.01, dynamic_binning_s=True, Single_tail_validation=True, tune_factor=0.99,
                     chisquaremethod=True, unimodal_categories=True):

        # GMM parameters
        tol = 1e-8  # convergence threshold
        n_mod = 3  # number of models for BIC comparison

        # BIC
        n_components = np.arange(1, n_mod)
        models = [GMM(n, tol=tol, random_state=41).fit(data)  # state=41 to make sure it doesn't always start random
                  for n in n_components]

        # BIC and AIC for validating optimal n
        BIC = [m.bic(data) for m in models]

        # use BIC as our standard due to limited data
        n_comp = np.argmin(BIC)+1

        # use guassian mixture modeling to model bimodal distribution
        gmm = GMM(n_components=n_comp, tol=tol, random_state=41).fit(data)

        # sort for cutoff
        x_axis2 = np.sort(np.array([b for a in data for b in a]))

        # recreate logscaled distribution
        # recreate normal distribution
        xs = np.linspace(np.min(data), np.max(data), num=499)
        xxs = np.linspace(np.min(data), np.max(data), num=500)
        dx = np.diff(xxs)  # calculate the integral value dx

        # recreate normal distribution with fitted GM
        means = gmm.means_
        covars = gmm.covariances_
        weights = gmm.weights_
        if n_comp > 1:
            yss = weights[0] * \
                stats.multivariate_normal.pdf(
                    xs, mean=means[0][0], cov=covars[0][0])
            yss2 = weights[1]*stats.multivariate_normal.pdf(
                xs, mean=means[1][0], cov=covars[1][0])
            expected = yss*dx*len(data)
            expectedt = yss2*dx*len(data)
            expectedboth = (yss+yss2)*dx*len(data)
            group_div = solve_gaussians(
                means[0][0], means[1][0], covars[0][0], covars[1][0], weights[0], weights[1], min(data), max(data))
        else:
            yss = stats.multivariate_normal.pdf(
                xs, mean=means[0][0], cov=covars[0][0])
            expected = yss*dx*len(data)
            expectedt = None
            expectedboth = None
            group_div = []

        # #finding out groups
        # groups = [np.argmax(a) for a in gmm.predict_proba(x_axis2.reshape(-1,1)).round(0)] #find which group each data point is under

        # #bimodality

        # roll = groups != np.roll(groups,1)
        # roll[0] = False  #roll will make the first variable True but we do not want that
        # group_index = [x for x in np.where(roll)][0]
        # group_div = x_axis2[group_index]

        # little distribution under big distribution (inclusion) is not the type of bimodality this algo is for return classif, categories, chis[-1]
        if len(group_div) == 0:
            n_comp = 1

    #     ------------------------------------measure chi-square of fitted GMM-------------------------------------------------
    #     determin number of bin using dynamic binning
        datas = np.array([[a] for a in np.sort(
            np.array([b for a in data for b in a]))])
        nums = int(1.88*len(datas)**(2/5))  # mann and wald

        observed, bin_edges = np.histogram(data, bins=np.linspace(
            np.min(data), np.max(data), num=nums), density=False)

        if dynamic_binning_s == True:
            golddof = len(observed)
            lenn = 0
            while lenn < golddof:  # Loop until you have the ssame dof or unimprovable

                observed, bin_edges, con = dynamic_binning(
                    observed, bin_edges, final=False)

                observed, what = np.histogram(data, bins=bin_edges)

                observed, bin_edges, con = dynamic_binning(observed, bin_edges)

                observed, what = np.histogram(data, bins=bin_edges)

                if con == 0:
                    break

                lenn = len(observed)

        # fit with a single GM
        # recreate normal distribution using dynamically binned edges
        xschi = bin_edges[1:]
        dxchi = np.diff(bin_edges)  # calculate the integral value dx

        if n_comp > 1:
            ysschi = weights[0]*stats.multivariate_normal.pdf(
                xschi, mean=means[0][0], cov=covars[0][0])
            yss2chi = weights[1]*stats.multivariate_normal.pdf(
                xschi, mean=means[1][0], cov=covars[1][0])
            expectedchi = (ysschi+yss2chi)*dxchi*len(data)
        else:
            ysschi = stats.multivariate_normal.pdf(
                xschi, mean=means[0][0], cov=covars[0][0])
            expectedchi = ysschi*dxchi*len(data)

        # calculate dof
        dof = len(observed)-1

        # calculate the critical value, *2 IS CHOSEN based on a priori and repeated testing
        critical_value = chi2.ppf(0.99999999, dof)

        # calculate chi-square
        arrchi = np.array([[x, y] for x, y in zip(observed, expectedchi)])

        cc = sum([(x[0]-x[1])**2./x[1] for x in arrchi])
        pp = 1-stats.chi2.cdf(cc, dof)

        bimodalc = cc

        bins = 10  # used for histogram and chi-square
        count = 0
        counting = 0
        conditions = [5]  # default with condition 5 meaning not fixed
        condition = 0  # default for those that doesn't go through any chi-square method
        chis = []
        rem = [0]  # default with condition 5 meaning not fixed so removes nothing
        BICs = []
        cclemon = []  # record chisquare
        biclemon = []  # record BIC
        realbic = []
        remove = None
        nums = int(1.88*len(data)**(2/5))  # mann and wald
        categories = None

    ########################################## UNIMODAL CATEGORIES ##########################################
    ########################################## UNIMODAL CATEGORIES ##########################################
    ########################################## UNIMODAL CATEGORIES ##########################################
    ########################################## UNIMODAL CATEGORIES ##########################################
    ########################################## UNIMODAL CATEGORIES ##########################################

        if unimodal_categories == True and n_comp == 1:
            # check if a small proportion of cells are background or expressing so we have two categories
            if verbosity == True:
                print('Unimodal distribution with distinct categories')

            twostd = meanf+2*stdf

            try:
                negpos_ratio = sum(data < filt)/len(data)
                posneg_ratio = sum(data > twostd)/len(data)
            except:
                raise Exception(
                    'Amigo, you did not put a background filter in, use the filt variable or make unimodal_categories = False if you do not need a background filter for the GMM pipeline')

            # if there is at least 5% of the cell lines that is in this category
            if (posneg_ratio <= 0.95 and posneg_ratio >= 0.05) and (negpos_ratio <= 0.95 and negpos_ratio >= 0.05):
                # set the threshold
                lower_bound = filt
                upper_bound = twostd

                categories = [1 if x <= lower_bound else 2 if x >=
                              upper_bound else np.nan for x in data]


    #     -----------------------------IF CHI-SQUARE VALUE > DOF WITH A SINGLE TAIL, RUN CHI-SQUARE METHOD-------
    #     -----------------------------IF CHI-SQUARE VALUE > DOF WITH A SINGLE TAIL, RUN CHI-SQUARE METHOD-------
    #     -----------------------------IF CHI-SQUARE VALUE > DOF WITH A SINGLE TAIL, RUN CHI-SQUARE METHOD-------
    #     -----------------------------IF CHI-SQUARE VALUE > DOF WITH A SINGLE TAIL, RUN CHI-SQUARE METHOD-------
    #     -----------------------------IF CHI-SQUARE VALUE > DOF WITH A SINGLE TAIL, RUN CHI-SQUARE METHOD-------
    #     -----------------------------IF CHI-SQUARE VALUE > DOF WITH A SINGLE TAIL, RUN CHI-SQUARE METHOD-------

        if Single_tail_validation == True:
            if cc > critical_value and n_comp == 2 and len(group_div) < 2:

                if verbosity == True:
                    print(
                        'Single tail problem that may not be normally distributed, run chi-square method')
                tail1 = [num for num in data if num <= group_div]
                tail2 = [num for num in data if num >= group_div]

                # deteremine which is the tail to begin trimming and fitting
                if len(tail1) > len(tail2):
                    thetail = tail2
                    direct = 0
                else:
                    thetail = tail1
                    direct = 1

                # determine the number of max datapoints to take off
                tailmax = math.ceil(len(thetail))
                if tailmax < 0.2*len(data):
                    tailmax = math.ceil(tailmax+farinto*tailmax)
                if verbosity == True:
                    print('Number of datapoints in tail:', tailmax)
                datas = np.array([[a] for a in np.sort(
                    np.array([b for a in data for b in a]))])  # sort data

                if len(thetail) > 0:  # if there are more than 1 datapoints in the tails then cut

                    while True:
                        for x in range(1, tailmax):

                            count = 0  # re-zero counter to 0 every loop
                            if direct == 0:
                                datamut = datas[:-x]
                            if direct == 1:
                                datamut = datas[x:]

                            nums = int(1.88*len(datamut)**(2/5))  # mann and wald

                            # BIC
                            models1 = [GMM(n, tol=tol, random_state=41).fit(datamut)
                                       for n in n_components]

                            # BIC and AIC for validating optimal n
                            BIC1 = [m.bic(datamut) for m in models1]

                            # use BIC as our standard due to limited data
                            n_comp1 = np.argmin(BIC1)+1

                            # dynamically making sure we result in the same number of bins
                            observed, bin_edges = np.histogram(datamut,
                                                               bins=np.linspace(np.min(datamut), np.max(datamut), num=nums), density=False)
                            if dynamic_binning_s == True:
                                golddof = len(observed)
                                lenn = 0
                                while lenn < golddof:

                                    observed, bin_edges, con = dynamic_binning(
                                        observed, bin_edges, final=False)

                                    observed, what = np.histogram(
                                        datamut, bins=bin_edges)

                                    observed, bin_edges, con = dynamic_binning(
                                        observed, bin_edges)

                                    observed, what = np.histogram(
                                        datamut, bins=bin_edges)

                                    if con == 0:
                                        break

                                    lenn = len(observed)

                            # fit with a single GM
                            gmmmut = GMM(n_components=1, tol=tol,
                                         random_state=41).fit(datamut)
                            # recreate normal distribution using dynamically binned edges
                            xs = bin_edges[1:]
                            # calculate the integral value dx
                            dx = np.diff(bin_edges)

                            # recreate normal distribution with fitted GM
                            means = gmmmut.means_
                            covars = gmmmut.covariances_
                            ys = stats.multivariate_normal.pdf(
                                xs, mean=means[0][0], cov=covars[0][0])
                            expected = dx*ys*np.sum(observed)

                            # calculate dof
                            dof = len(observed)-1

                            # calculate chi-square
                            arr = np.array([[x, y]
                                            for x, y in zip(observed, expected)])
                            c = sum([(x[0]-x[1])**2./x[1] for x in arr])
                            p = 1-stats.chi2.cdf(c, dof)

                            # fit with two GMs
                            gmmmut2 = GMM(n_components=2, tol=tol,
                                          random_state=41).fit(datamut)

                            # figure where there are two groups still or not
                            g1 = [np.argmax(a) for a in gmmmut2.predict_proba(
                                datamut.reshape(-1, 1)).round(0)]  # find which group each data point is under
                            # bimodality
                            r1 = g1 != np.roll(g1, 1)
                            # roll will make the first variable True but we do not want that
                            r1[0] = False
                            gi1 = datamut[[x for x in np.where(r1)]]

                            # recreate normal distribution with 2 fitted GM
                            weights2 = gmmmut2.weights_
                            means2 = gmmmut2.means_
                            covars2 = gmmmut2.covariances_

                            yss = weights2[0]*stats.multivariate_normal.pdf(
                                xs, mean=means2[0][0], cov=covars2[0][0])
                            yss2 = weights2[1]*stats.multivariate_normal.pdf(
                                xs, mean=means2[1][0], cov=covars2[1][0])
                            expected2 = (yss+yss2)*dx*np.sum(observed)

                            # calculate chi-square
                            arr2 = np.array([[x, y]
                                             for x, y in zip(observed, expected2)])
                            c2 = sum([(x[0]-x[1])**2./x[1] for x in arr2])
                            p2 = 1-stats.chi2.cdf(c2, dof)

                            # reset xs
                            # recreate normal distribution
                            xs = np.linspace(np.min(data), np.max(data), num=499)
                            xxs = np.linspace(np.min(data), np.max(data), num=500)
                            dx = np.diff(xxs)  # calculate the integral value dx

                            # is it better than the original fit?
                            if counting == 0:
                                # degrees of freedom factor
                                ctf = round(cc*tune_factor, 2)
                                if x == 1:
                                    cclemon.append([ctf, ctf])
                                    biclemon.append([n_comp])
                                    realbic.append(BIC1)
                                    chis.append(ctf)
                                    rem.append(1)
                            else:
                                # chisquare tunning factor
                                ctf = round(chis[-1]*tune_factor, 2)
                            if n_comp1 == 1:
                                if c < ctf:
                                    if verbosity == True:
                                        print(
                                            'Removed %d datanormpoints and fitting with 1 Gaussian Mixture with Chi-square value %s <= %s' % (x, np.round(c, 2), ctf))
                                    count = 1
                                    fc = np.round(c, 2)

                            if n_comp1 == 2:
                                if len(gi1) < 2:
                                    if c2 <= ctf/2:
                                        if verbosity == True:
                                            print('Removed %d datanormpoints and fitting with 2 Gaussian Mixture with Chi-square value %s <= %s' % (
                                                x, np.round(c2, 2), ctf))
                                        count = 2
                                        fc = np.round(c2, 2)

                            if count > 0:  # only begin recording when it fits under one of the conditions above
                                if direct == 1:
                                    conditions.append(count)
                                    chis.append(fc)
                                    counting += 1
                                    rem.append(x)
                                    BICs.append(BIC1)

                                if direct == 0:
                                    # if direct == 0 means we are trimming from the right so add 2 for conditions and count
                                    conditions.append(count+2)
                                    chis.append(fc)
                                    counting += 1
                                    rem.append(x)
                                    BICs.append(BIC1)
                                    if n_comp1 == 1:
                                        count = 3
                                    if n_comp1 == 2:
                                        count = 4

                            cclemon.append([c, c2])
                            biclemon.append([n_comp1])
                            realbic.append(BIC1)

                        if count == 0:
                            # fifth condition where it did not fix
                            conditions = [5]
                            break

                    condition = conditions[-1]
                    remove = rem[-1]

                    if condition == 1:
                        datamut = datas[remove:]
                        n_comp = 1
                        gmm = GMM(n_components=1, tol=tol,
                                  random_state=41).fit(datamut)
                    elif condition == 2:
                        datamut = datas[remove:]
                        n_comp = 2
                        gmm = GMM(n_components=2, tol=tol,
                                  random_state=41).fit(datamut)
                    elif condition == 3:
                        datamut = datas[:-remove]
                        n_comp = 1
                        gmm = GMM(n_components=1, tol=tol,
                                  random_state=41).fit(datamut)
                    elif condition == 4:
                        datamut = datas[:-remove]
                        n_comp = 2
                        gmm = GMM(n_components=2, tol=tol,
                                  random_state=41).fit(datamut)
                    elif condition == 5:
                        n_comp = 2
                        gmm = GMM(n_components=2, tol=tol,
                                  random_state=41).fit(data)
                        if verbosity == True:
                            print('Chi-square Method Did Not Help the Tail Problem')

                    # BIC
                    if condition == 1 or condition == 2 or condition == 3 or condition == 4:
                        BIC = BICs[-1]
                    else:
                        n_components = np.arange(1, n_mod)
                        models = [GMM(n, tol=tol, random_state=41).fit(data)
                                  for n in n_components]

                        # BIC and AIC for validating optimal n
                        BIC = [m.bic(data) for m in models]

                        # optimal n_comp
                        n_comp = np.argmin(BIC)+1

                    # use guassian mixture modeling to model bimodal distribution
                    if condition == 1 or condition == 2 or condition == 3 or condition == 4:
                        x_axis2 = np.sort(
                            np.array([b for a in datamut for b in a]))
                    else:
                        x_axis2 = np.sort(np.array([b for a in data for b in a]))

                    # recreate normal distribution with fitted GM
                    means = gmm.means_
                    covars = gmm.covariances_
                    weights = gmm.weights_

                    if n_comp > 1:
                        yss = weights[0]*stats.multivariate_normal.pdf(
                            xs, mean=means[0][0], cov=covars[0][0])
                        yss2 = weights[1]*stats.multivariate_normal.pdf(
                            xs, mean=means[1][0], cov=covars[1][0])
                        expected = yss*dx*len(data)
                        expectedt = yss2*dx*len(data)
                        expectedboth = (yss+yss2)*dx*len(data)
                        group_div = solve_gaussians(
                            means[0][0], means[1][0], covars[0][0], covars[1][0], weights[0], weights[1], min(data), max(data))
                    else:
                        yss = stats.multivariate_normal.pdf(
                            xs, mean=means[0][0], cov=covars[0][0])
                        expected = yss*dx*len(data)
                        expectedt = None
                        expectedboth = None
                        group_div = []

                    # #finding out groups
                    # groups = [np.argmax(a) for a in gmm.predict_proba(x_axis2.reshape(-1,1)).round(0)] #find which group each data point is under

                    # #bimodality

                    # roll = groups != np.roll(groups,1)
                    # roll[0] = False  #roll will make the first variable True but we do not want that
                    # group_index = [x for x in np.where(roll)][0]
                    # group_div = x_axis2[group_index]
                else:
                    n_comp = 1

    #     ----------------------------IF TAIL PROBELM RUN ------------------------------------------------------------
        if chisquaremethod == True:
            if len(group_div) > 1:
                if verbosity == True:
                    print('Rerunning GMM with Chi-square Method to fix tail problem')

                tail1 = [num for num in data if num <= np.min(group_div)]
                tail2 = [num for num in data if num >= np.max(group_div)]

                chiv = []
                pv = []
                chiv2 = []
                pv2 = []
                chiv3 = []
                pv3 = []
                chiv4 = []
                pv4 = []
                BICs = []
                gs = []
                gs2 = []

                # determine the number of max datapoints to take off
                xup = math.ceil(len(tail1))
                yup = math.ceil(len(tail2))

                # add to the boundary if the tail is very small to allow more data taken off for validation
                if xup < 0.2*len(data):
                    xup = math.ceil(xup+xup*farinto)
                if yup < 0.2*len(data):
                    # plus two to avoid 0 len tails
                    yup = math.ceil(yup+yup*farinto)

                if verbosity == True:
                    print('Number of datapoints in the right tail:', xup)
                    print('Number of datapoints in the left tail:', yup)

                while True:
                    for x, y in itert.zip_longest(range(1, xup), range(1, yup)):

                        # sort data so we are taking off the right tail
                        datas = np.array([[a] for a in np.sort(
                            np.array([b for a in data for b in a]))])
                        count = 0  # re-zero counter to 0 every loop
                        try:
                            datamut = datas[x:]  # tail1 data
                            datamut2 = datas[:-y]  # tail2 data

                        except:
                            pass

                        nums = int(1.88*len(datamut)**(2/5))  # mann and wald

                        # BIC
                        models1 = [GMM(n, tol=tol, random_state=41).fit(datamut)
                                   for n in n_components]

                        # BIC and AIC for validating optimal n
                        BIC1 = [m.bic(datamut) for m in models1]

                        # use BIC as our standard due to limited data
                        n_comp1 = np.argmin(BIC1)+1

                    # ---------------------------------tail1----------------------------------------------------

                        observed, bin_edges = np.histogram(datamut,
                                                           bins=np.linspace(np.min(datamut), np.max(datamut), num=nums), density=False)
                        if dynamic_binning_s == True:
                            golddof = len(observed)
                            lenn = 0
                            while lenn < golddof:
                                observed, bin_edges, con = dynamic_binning(
                                    observed, bin_edges, final=False)

                                observed, what = np.histogram(
                                    datamut, bins=bin_edges)

                                observed, bin_edges, con = dynamic_binning(
                                    observed, bin_edges)

                                observed, what = np.histogram(
                                    datamut, bins=bin_edges)

                                if con == 0:
                                    break

                                lenn = len(observed)
                        # fit with a single GM
                        gmmmut = GMM(n_components=1, tol=tol,
                                     random_state=41).fit(datamut)
                        # recreate normal distribution using dynamically binned edges
                        xs = bin_edges[1:]
                        dx = np.diff(bin_edges)  # calculate the integral value dx

                        # recreate normal distribution with fitted GM
                        means = gmmmut.means_
                        covars = gmmmut.covariances_
                        ys = stats.multivariate_normal.pdf(
                            xs, mean=means[0][0], cov=covars[0][0])
                        expected = dx*ys*np.sum(observed)

                        # calculate dof
                        dof = len(observed)-1

                        # calculate chi-square
                        arr = np.array([[x, y]
                                        for x, y in zip(observed, expected)])
                        c = sum([(x[0]-x[1])**2./x[1] for x in arr])
                        p = 1-stats.chi2.cdf(c, dof)
                        chiv.append(c), pv.append(p)

                        # fit with two GMs
                        gmmmut2 = GMM(n_components=2, tol=tol,
                                      random_state=41).fit(datamut)

                        # figure where there are two groups still or not
                        g1 = [np.argmax(a) for a in gmmmut2.predict_proba(
                            datamut.reshape(-1, 1)).round(0)]  # find which group each data point is under
                        # bimodality
                        r1 = g1 != np.roll(g1, 1)
                        # roll will make the first variable True but we do not want that
                        r1[0] = False
                        gi1 = datamut[[x for x in np.where(r1)]]
                        gs.append(len(gi1))

                        # recreate normal distribution with 2 fitted GM
                        weights2 = gmmmut2.weights_
                        means2 = gmmmut2.means_
                        covars2 = gmmmut2.covariances_
                        yss = weights2[0]*stats.multivariate_normal.pdf(
                            xs, mean=means2[0][0], cov=covars2[0][0])
                        yss2 = weights2[1]*stats.multivariate_normal.pdf(
                            xs, mean=means2[1][0], cov=covars2[1][0])
                        expected2 = (yss+yss2)*dx*np.sum(observed)

                        # calculate chi-square
                        arr2 = np.array([[x, y]
                                         for x, y in zip(observed, expected2)])
                        c2 = sum([(x[0]-x[1])**2./x[1] for x in arr2])
                        p2 = 1-stats.chi2.cdf(c2, dof)
                        chiv2.append(c2), pv2.append(p2)

                        # BIC
                        models2 = [GMM(n, tol=tol, random_state=41).fit(datamut2)
                                   for n in n_components]

                        # BIC and AIC for validating optimal n
                        BIC2 = [m.bic(datamut2) for m in models2]

                        # use BIC as our standard due to limited data
                        n_comp2 = np.argmin(BIC2)+1
                    # ---------------------------------tail2----------------------------------------------------
                        observed2, bin_edges2 = np.histogram(datamut2,
                                                             bins=np.linspace(np.min(datamut2), np.max(datamut2), num=nums), density=False)

                        if dynamic_binning_s == True:
                            golddof = len(observed2)
                            lenn = 0
                            while len(observed2) < golddof:
                                observed2, bin_edges2, con = dynamic_binning(
                                    observed2, bin_edges2, final=False)

                                observed2, what = np.histogram(
                                    datamut2, bins=bin_edges2)

                                observed2, bin_edges2, con = dynamic_binning(
                                    observed2, bin_edges2)

                                observed2, what = np.histogram(
                                    datamut2, bins=bin_edges2)

                                if con == 0:
                                    lenn = len(observed2)

                        # fit with a single GM
                        gmmmut3 = GMM(n_components=1, tol=tol,
                                      random_state=41).fit(datamut2)
                        xs2 = bin_edges2[1:]  # recreate normal distribution
                        # calculate the integral value dx
                        dx2 = np.diff(bin_edges2)

                        # recreate normal distribution with fitted GM cutting from tail2
                        meanso = gmmmut3.means_
                        covarso = gmmmut3.covariances_
                        yso = stats.multivariate_normal.pdf(
                            xs2, mean=meanso[0][0], cov=covarso[0][0])
                        expectedo = dx2*yso*np.sum(observed2)

                        # calculate chi-square
                        arro = np.array([[x, y]
                                         for x, y in zip(observed2, expectedo)])
                        c3 = sum([(x[0]-x[1])**2./x[1] for x in arro])
                        p3 = 1-stats.chi2.cdf(c3, dof)
                        chiv3.append(c3), pv3.append(p3)
                        # fit with two GMs
                        gmmmut4 = GMM(n_components=2, tol=tol,
                                      random_state=41).fit(datamut2)

                        # figure where there are two groups still or not
                        g2 = [np.argmax(a) for a in gmmmut4.predict_proba(
                            datamut2.reshape(-1, 1)).round(0)]  # find which group each data point is under
                        # bimodality
                        r2 = g2 != np.roll(g2, 1)
                        # roll will make the first variable True but we do not want that
                        r2[0] = False
                        gi2 = datamut2[[x for x in np.where(r2)]]
                        gs2.append(len(gi2))

                        # recreate normal distribution with 2 fitted GM from tail2
                        weightso2 = gmmmut4.weights_
                        meanso2 = gmmmut4.means_
                        covarso2 = gmmmut4.covariances_
                        ysso = weightso2[0]*stats.multivariate_normal.pdf(
                            xs2, mean=meanso2[0][0], cov=covarso2[0][0])
                        ysso2 = weightso2[1]*stats.multivariate_normal.pdf(
                            xs2, mean=meanso2[1][0], cov=covarso2[1][0])
                        expectedo2 = (ysso+ysso2)*dx2*np.sum(observed2)

                        # calculate chi-square
                        arro2 = np.array([[x, y]
                                          for x, y in zip(observed2, expectedo2)])
                        c4 = sum([(x[0]-x[1])**2./x[1] for x in arro2])
                        p4 = 1-stats.chi2.cdf(c4, dof)
                        chiv4.append(c4), pv4.append(p4)

                        # reset xs
                        # recreate normal distribution
                        xs = np.linspace(np.min(data), np.max(data), num=499)
                        xxs = np.linspace(np.min(data), np.max(data), num=500)
                        dx = np.diff(xxs)  # calculate the integral value dx

                        if counting == 0:
                            # degrees of freedom factor
                            ctf = round(cc*tune_factor, 2)
                            if x == 1:
                                chis.append(ctf)
                                cclemon.append([ctf, ctf, ctf, ctf])
                                biclemon.append([n_comp, n_comp])
                                realbic.append([BIC1, BIC2])
                                rem.append(1)
                        else:
                            # chisquare tunning factor
                            ctf = round(chis[-1]*tune_factor, 2)
                        # stop when p value is lower than <0.05 , find what condition?
                        fc = 0

                        if n_comp1 == 1:
                            if x != None:
                                if c < ctf:
                                    if verbosity == True:
                                        print(
                                            'Removed %d datapoints from the left tail and fitting with 1 Gaussian Mixture with Chi-square value %s <= %s' % (x, np.round(c, 2), ctf))
                                    count = 1
                                    fc = np.round(c, 2)
                        if n_comp2 == 1:
                            if y != None:
                                if c3 < ctf:
                                    if verbosity == True:
                                        print('Removed %d datapoints from the right tail and fitting with 1 Gaussian Mixture with Chi-square value %s <= %s' % (
                                            y, np.round(c3, 2), ctf))
                                    count = 3
                                    fc = np.round(c3, 2)
                        if n_comp1 == 2:
                            if x != None:
                                if len(gi1) < 2:
                                    if c2 < ctf/2:
                                        if verbosity == True:
                                            print('Removed %d datapoints from the left tail and fitting with 2 Gaussian Mixture with Chi-square value %s <= %s' % (
                                                x, np.round(c2, 2), ctf))
                                        count = 2
                                        fc = np.round(c2, 2)
                        if n_comp2 == 2:
                            if y != None:
                                if len(gi2) < 2:
                                    if c4 < ctf/2:
                                        if verbosity == True:
                                            print('Removed %d datapoints from the right tail and fitting with 2 Gaussian Mixture with Chi-square value %s <= %s' % (
                                                y, np.round(c4, 2), ctf))
                                        count = 4
                                        fc = np.round(c4, 2)

                        cclemon.append([c, c2, c3, c4])
                        biclemon.append([n_comp1, n_comp2])
                        realbic.append([BIC1, BIC2])

                        if count > 0:  # only begin recording when it fits under one of the conditions above
                            conditions.append(count)
                            chis.append(fc)
                            counting += 1
                        if count == 1 or count == 2:
                            rem.append(x)
                            BICs.append(BIC1)
                        elif count == 3 or count == 4:
                            rem.append(y)
                            BICs.append(BIC2)

                    if count == 0:
                        conditions = [5]  # fifth condition where it did not fix
                        break

                # -----------------------------------------------rerun GMM----------------------------------------------------
                condition = conditions[-1]
                remove = rem[-1]

                if condition == 1:
                    datamut = datas[remove:]
                    n_comp = 1
                    gmm = GMM(n_components=1, tol=tol,
                              random_state=41).fit(datamut)
                elif condition == 2:
                    datamut = datas[remove:]
                    n_comp = 2
                    gmm = GMM(n_components=2, tol=tol,
                              random_state=41).fit(datamut)
                elif condition == 3:
                    datamut = datas[:-remove]
                    n_comp = 1
                    gmm = GMM(n_components=1, tol=tol,
                              random_state=41).fit(datamut)
                elif condition == 4:
                    datamut = datas[:-remove]
                    n_comp = 2
                    gmm = GMM(n_components=2, tol=tol,
                              random_state=41).fit(datamut)
                elif condition == 5:
                    n_comp = 2
                    gmm = GMM(n_components=2, tol=tol, random_state=41).fit(data)
                    if verbosity == True:
                        print('Chi-square Method Did Not Fix the Tail Problem')

                # BIC
                if condition == 1 or condition == 2 or condition == 3 or condition == 4:
                    BIC = BICs[-1]
                else:
                    n_components = np.arange(1, n_mod)
                    models = [GMM(n, tol=tol, random_state=41).fit(data)
                              for n in n_components]

                    # BIC and AIC for validating optimal n
                    BIC = [m.bic(data) for m in models]

                    # optimal n_comp
                    n_comp = np.argmin(BIC)+1

                # use guassian mixture modeling to model bimodal distribution
                if condition == 1 or condition == 2:
                    x_axis2 = np.sort(np.array([b for a in datamut for b in a]))
                elif condition == 3 or condition == 4:
                    x_axis2 = np.sort(np.array([b for a in datamut for b in a]))
                else:
                    x_axis2 = np.sort(np.array([b for a in data for b in a]))

                # recreate normal distribution with fitted GM
                means = gmm.means_
                covars = gmm.covariances_
                weights = gmm.weights_
                if n_comp > 1:
                    yss = weights[0]*stats.multivariate_normal.pdf(
                        xs, mean=means[0][0], cov=covars[0][0])
                    yss2 = weights[1]*stats.multivariate_normal.pdf(
                        xs, mean=means[1][0], cov=covars[1][0])
                    expected = yss*dx*len(data)
                    expectedt = yss2*dx*len(data)
                    expectedboth = (yss+yss2)*dx*len(data)
                    group_div = solve_gaussians(
                        means[0][0], means[1][0], covars[0][0], covars[1][0], weights[0], weights[1], min(data), max(data))
                else:
                    yss = stats.multivariate_normal.pdf(
                        xs, mean=means[0][0], cov=covars[0][0])
                    expected = yss*dx*len(data)
                    expectedt = None
                    expectedboth = None
                    group_div = []

                # #finding out groups
                # groups = [np.argmax(a) for a in gmm.predict_proba(x_axis2.reshape(-1,1)).round(0)] #find which group each data point is under

                # #bimodality

                # roll = groups != np.roll(groups,1)
                # roll[0] = False  #roll will make the first variable True but we do not want that
                # group_index = [x for x in np.where(roll)][0]
                # group_div = x_axis2[group_index]

        return weights, categories, group_div, datas, condition, chis, rem, BICs, cclemon, biclemon, realbic, n_comp, n_components, BIC, bimodalc, xs, expected, expectedt, expectedboth, nums, count, remove, dof, means, covars


    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    # v
    ##############################################################################################################################
    ##############################################################################################################################


    def GMM_plot(input_data, data, datanorm, input_datanorm, categories, verbosity, graphs, graph_show, n_comp, xs, expected, expectedt, expectedboth, ID, log2transform, nums, meanf,
                 stdf, group_div, count, condition, n_components, BIC, calc_back, input_datanormcat,
                 chis, means, covars, bimodalc, remove, datas, biclemon, realbic, filt, cell_lines, cclemon, dof,
                 rem, cell_line_groupname, output_bins):
        from matplotlib.ticker import MaxNLocator

       # -------------------------------------------PLOTTTTTTTTTTTTTTTTT------------------------------------------
        x_axis = np.sort(np.array([b for a in datanorm[np.nonzero(
            datanorm)].reshape(-1, 1) for b in a]))  # reset x_axis for plotting
        # initialize f
        f = 0
        if graphs == True:
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=[18, 14.5])
            ax1.yaxis.set_major_locator(MaxNLocator(integer=True))

            axisfont = 18
            titlefont = 18.5

            a2 = ax1.twinx()
            # plot
            if log2transform == True:  # 2** if log2 transformed to go it back
                xs = 2**xs
            elif log2transform == False:
                xs = xs

            if n_comp > 1:
                a2.plot(xs, expected, '--k', color='black')
                a2.plot(xs, expectedt, '--k', color='black')
                a2.plot(xs, expectedboth, '-k', color='b')
            else:
                a2.plot(xs, expected, '-k', color='b')

            a2.set_ylim(bottom=0)
            a2.set_yticks([])

            if log2transform == True:

                # plot histogram of original data
                ax1.set_title(ID, fontsize=titlefont)

                n, bins, patches = ax1.hist(datanorm,
                                            bins=np.logspace(np.log2(np.min(x_axis)), np.log2(np.max(x_axis)), num=nums, base=2), histtype='barstacked', alpha=0.5, ec='black')
                ax1.set_ylabel('Number of Samples', fontsize=axisfont)
                ax1.set_xlabel('Expression Value (log2)', fontsize=axisfont)
                ax1.set_xscale('log', basex=2)

                # take out the 2^ just have the powers
                logfmt = tick.LogFormatterExponent(base=2.0, labelOnlyBase=True)
                ax1.xaxis.set_major_formatter(logfmt)
                ax1.tick_params(axis='both', labelsize=axisfont)
                a2.tick_params(axis='both', labelsize=axisfont)
                ax1.grid('on')

                if filt != None:
                    twostd = 2**(meanf+2*stdf)
                    filt = 2**filt
                else:
                    filt = -10000000000
                    # make it the smallest data so that it doesn't interefere with graphing
                    twostd = -10000000000

            elif log2transform == False:
                # plot histogram of original data
                ax1.set_title(ID, fontsize=titlefont)

                n, bins, patches = ax1.hist(
                    datanorm, nums, histtype='barstacked', alpha=0.5, ec='black')
                ax1.set_ylabel('Number of Samples', fontsize=axisfont)
                ax1.set_xlabel('Expression Value', fontsize=axisfont)

                ax1.tick_params(axis='both', labelsize=axisfont)
                a2.tick_params(axis='both', labelsize=axisfont)
                ax1.grid('on')

                if filt != None:
                    twostd = meanf+2*stdf
                else:
                    filt = -10000000000
                    # make it the smallest data so that it doesn't interefere with graphing
                    twostd = -10000000000

            # color in background threshold
            if calc_back == False and filt != None:

                # color potential background, x < threshold
                if sum(datanorm < twostd) > 1:
                    ind = [i for i, y in enumerate(bins) if y <= twostd]
                    if len(ind) > 1:
                        for i in range(ind[0], ind[-1]):
                            patches[i].set_facecolor('orange')
                if np.min(datanorm) < filt:
                    ind = [i for i, y in enumerate(bins) if y <= filt]
                    if len(ind) > 1:
                        for i in range(ind[0], ind[-1]):
                            patches[i].set_facecolor('r')

                background_p = mpatches.Patch(
                    facecolor='lightcoral', ec='black', label='Background Threshold')
                background_p2 = mpatches.Patch(
                    facecolor='gold', ec='black', label='Two \u03C3 Background Threshold')

                a2.legend(handles=[background_p, background_p2])

            # ADD PATCHES AND HATCHES FOR USER INPUT CELL_LINES
            if cell_lines != []:
                cmap = plt.get_cmap('tab20')  # color scheme

                all_highlight = []

                # check if list of list or just a list of cell lines
                if all([isinstance(i, list) for i in cell_lines]):
                    # stack the hatches
                    stacked = []
                    stacked_color = []

                    for i, x in enumerate(cell_lines):

                        use_color = cmap(i/len(cell_lines))  # cycle through color

                        highlight = input_datanorm[x].loc[ID]
                        if cell_line_groupname != []:  # add groupname into the patches to show
                            try:
                                all_highlight.append(mpatches.Patch(hatch='//', facecolor=use_color, edgecolor='black',
                                                                    label=cell_line_groupname[i]))
                            except:
                                print(
                                    'cell_line_groupname needs to be the same length as cell_lines')
                                all_highlight.append(mpatches.Patch(hatch='//', facecolor=use_color, edgecolor='black',
                                                                    label='Sample(s) of Interest No.{}'.format(i)))  # add to legend
                        else:
                            all_highlight.append(mpatches.Patch(hatch='//', facecolor=use_color, edgecolor='black',
                                                                label='Sample(s) of Interest No.{}'.format(i)))  # add to legend

                        stacked.append(np.array(input_datanorm[x].loc[ID].values))
                        stacked_color.append(use_color)

                    # bin spaces different for log2transformed or not
                    if log2transform == True:
                        n2, bins2, patches2 = ax1.hist(stacked, hatch='//', color=stacked_color,
                                                       bins=np.logspace(np.log2(np.min(x_axis)), np.log2(
                                                           np.max(x_axis)), num=nums, base=2),
                                                       stacked=True, alpha=0.7, ec='black')
                    else:
                        n2, bins2, patches2 = ax1.hist(stacked, hatch='//', color=stacked_color,
                                                       bins=np.linspace(min(x_axis), max(x_axis), num=nums+1), stacked=True,
                                                       alpha=0.7, ec='black')

                else:  # if just a list of cell lines then just hatch everything as if one group
                    highlight = input_datanorm[cell_lines].loc[ID]

                    all_highlight = [mpatches.Patch(hatch='///', facecolor='w', edgecolor='black',
                                                    label='Sample(s) of Interest')]  # add to legend

                    if log2transform == True:
                        n2, bins2, patches2 = ax1.hist(input_datanorm[cell_lines].loc[ID], hatch='///', facecolor='none',
                                                       bins=np.logspace(np.log2(np.min(x_axis)), np.log2(
                                                           np.max(x_axis)), num=nums, base=2),
                                                       histtype='barstacked', alpha=0.5, ec='black')
                    else:
                        n2, bins2, patches2 = ax1.hist(input_datanorm[cell_lines].loc[ID], hatch='///', facecolor='none',
                                                       bins=np.linspace(
                                                           min(x_axis), max(x_axis), num=nums+1),
                                                       histtype='barstacked', alpha=0.5, ec='black')

                a2.legend(handles=[background_p, background_p2] +
                          all_highlight, loc='upper right')  # add all labels

            if log2transform == True:
                # deem the output unimodal or bimodal whether chi-square method was applied
                intersections = [2**x for x in group_div]

                if count == 0:
                    if n_comp > 1:
                        for a in intersections:
                            ax1.axvline(a, linestyle='--', c='red')
                        if verbosity == True:
                            print('Bimodal, Cutoff Threshold:', group_div)
                        xf = group_div
                        classif = 'bimodal'
                    else:
                        pass
                        if verbosity == True:
                            print('Unimodal')
                        xf = None
                        if categories != None:
                            classif = 'categorical unimodal'
                        else:
                            classif = 'unimodal'
                elif condition == 1:  # remove 1 point means removeing 0 index need to zero index
                    ax1.axvline(2**datas[remove][0], linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the left tail, Two Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, datas[remove][0]))
                    xf = datas[remove][0]
                    classif = 'unimodal with a tail'
                elif condition == 2:
                    ax1.axvline(2**datas[remove][0], linestyle='--', c='red')
                    for a in intersections:
                        ax1.axvline(a, linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the left tail, Three Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, np.sort(np.concatenate([group_div, datas[remove]]))))
                    xf = np.sort(np.concatenate([group_div, datas[remove]]))
                    classif = 'bimodal with a tail'
                elif condition == 3:
                    remd = remove+1  # not zero indexed coming from the other side
                    ax1.axvline(2**datas[-remd][0], linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the right tail, Two Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, datas[-remove][0]))
                    xf = datas[-remove][0]
                    classif = 'unimodal with a tail'
                elif condition == 4:
                    remd = remove+1
                    ax1.axvline(2**datas[-remd][0], linestyle='--', c='red')
                    for a in intersections:
                        ax1.axvline(a, linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the right tail, Three Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, np.sort(np.concatenate([group_div, datas[-remove]]))))
                    xf = np.sort(np.concatenate([group_div, datas[-remove]]))
                    classif = 'bimodal with a tail'
                elif condition == 5:
                    for a in intersections:
                        ax1.axvline(a, linestyle='--', c='red')
                    if verbosity == True:
                        print(
                            'Bimodal cannot be helped by Chi-square Method, Cutoff Threshold:', group_div)
                    xf = group_div
                    classif = 'poorly fitted bimodal'

            elif log2transform == False:
                # deem the output unimodal or bimodal whether chi-square method was applied
                intersections = group_div

                if count == 0:
                    if n_comp > 1:
                        for a in intersections:
                            ax1.axvline(a, linestyle='--', c='red')
                        if verbosity == True:
                            print('Bimodal, Cutoff Threshold:', group_div)
                        xf = group_div
                        classif = 'bimodal'
                    else:
                        pass
                        if verbosity == True:
                            print('Unimodal')
                        xf = None
                        if categories != None:
                            classif = 'categorical unimodal'
                        else:
                            classif = 'unimodal'
                elif condition == 1:  # remove 1 point means removeing 0 index need to zero index
                    ax1.axvline(datas[remove][0], linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the left tail, Two Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, datas[remove][0]))
                    xf = datas[remove][0]
                    classif = 'unimodal with a tail'
                elif condition == 2:
                    ax1.axvline(datas[remove][0], linestyle='--', c='red')
                    for a in intersections:
                        ax1.axvline(a, linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the left tail, Three Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, np.sort(np.concatenate([group_div, datas[remove]]))))
                    xf = np.sort(np.concatenate([group_div, datas[remove]]))
                    classif = 'bimodal with a tail'
                elif condition == 3:
                    remd = remove+1  # not zero indexed coming from the other side
                    ax1.axvline(datas[-remd][0], linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the right tail, Two Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, datas[-remove][0]))
                    xf = datas[-remove][0]
                    classif = 'unimodal with a tail'
                elif condition == 4:
                    remd = remove+1
                    ax1.axvline(datas[-remd][0], linestyle='--', c='red')
                    for a in intersections:
                        ax1.axvline(a, linestyle='--', c='red')
                    if verbosity == True:
                        print('Removed %s datapoints from the right tail, Three Categories by Chi-square Method, Cutoff Threshold: %s' %
                              (remove, np.sort(np.concatenate([group_div, datas[-remove]]))))
                    xf = np.sort(np.concatenate([group_div, datas[-remove]]))
                    classif = 'bimodal with a tail'
                elif condition == 5:
                    for a in intersections:
                        ax1.axvline(a, linestyle='--', c='red')
                    if verbosity == True:
                        print(
                            'Bimodal cannot be helped by Chi-square Method, Cutoff Threshold:', group_div)
                    xf = group_div
                    classif = 'poorly fitted bimodal'

        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        # BIC and AIC graph
            ax2.set_title('BIC Graph of :' + ' ' + ID, fontsize=titlefont)
            ax2.tick_params(axis='both', labelsize=axisfont)

            if condition != 0 and condition != 5:
                xp = np.arange(0, len(biclemon))
                length = len(biclemon[0])
                ax2.set_ylim((0.8, 2.2))
                ax2.set_yticks([1, 2])
                ax2.set_ylabel('n_components', fontsize=axisfont)
                ax2.set_xlabel('Datapoints Removed', fontsize=axisfont)

                ax22 = ax2.twinx()  # plot BIC value
                ax22.tick_params(axis='both', labelsize=axisfont)

                if length == 1:
                    ax2.plot(xp, biclemon, '--', c='grey')

                    color = ['b', 'r']
                    for zz in range(0, len(realbic[0])):
                        bicbag = []
                        for bics in realbic:
                            bicbag.append(bics[zz])
                        ax22.plot(xp, bicbag, '-', c=color[zz])

                    b1 = mlines.Line2D([], [], color=color[0],
                                       label='1 n_comp w/ BIC')
                    b2 = mlines.Line2D([], [], color=color[1],
                                       label='2 n_comp w/ BIC')
                    b3 = mlines.Line2D([], [], color='grey',
                                       label='n_comp w/ Lower BIC Score')
                    ax22.legend(handles=[b1, b2, b3])
                    ax22.set_ylabel('BIC Score', fontsize=axisfont)
                else:
                    color = ['grey', 'black']
                    for zz in range(0, length):
                        bicbags = []
                        for bics in biclemon:
                            bicbags.append(bics[zz])
                        ax2.plot(xp, bicbags, '--', c=color[zz])

                    colors = ['deepskyblue', 'dodgerblue', 'tomato', 'lightcoral']

                    newbic = [np.concatenate(x) for x in realbic]
                    for zz in range(0, len(newbic[0])):
                        bicbag = []
                        for bics in newbic:
                            bicbag.append(bics[zz])
                        ax22.plot(xp, bicbag, '-', c=colors[zz])

                    bb1 = mlines.Line2D([], [], color=colors[0],
                                        label='Left Tail w/ 1 n_comp BIC')
                    bb2 = mlines.Line2D([], [], color=colors[1],
                                        label='Left Tail w/ 2 n_comp BIC')
                    bb3 = mlines.Line2D([], [], color=colors[2],
                                        label='Right Tail w/ 1 n_comp BIC')
                    bb4 = mlines.Line2D([], [], color=colors[3],
                                        label='Right Tail w/ 2 n_comp BIC')
                    b1 = mlines.Line2D([], [], color=color[0],
                                       label='Left tail n_comp w/ Lower BIC')
                    b2 = mlines.Line2D([], [], color=color[1],
                                       label='Right tail n_comp w/ Lower BIC')
                    ax22.legend(handles=[bb1, bb2, bb3, bb4, b1, b2])
                    ax22.set_ylabel('BIC Score', fontsize=axisfont)
            else:
                ax2.bar(n_components, BIC, width=0.2,
                        label='BIC', align='center', ec='b')
                ax2.set_xticks([1, 2])
                ax2.set_ylim(bottom=np.min(BIC)-5, top=np.max(BIC)+2)
                ax2.set_ylabel('BIC Score', fontsize=axisfont)
                ax2.set_xlabel('n_components', fontsize=axisfont)
            # ------------------------------------------------------------------------------------------------------
            # ------------------------------------------------------------------------------------------------------
            # ------------------------------------------------------------------------------------------------------
            # ------------------------------------------------------------------------------------------------------

            # plot number of samples in each category

            ax3.grid('on')
            ax3.tick_params(axis='both', labelsize=axisfont)
            cutoff = xf
            try:
                if len(cutoff) == 1:
                    cutoff = cutoff[0]
            except:
                pass

            # reset twostd
            twostd = meanf+2*stdf

            if sum(cutoff == None) == 1:
                if categories == None:
                    categories = np.ones(input_data.shape[1]).tolist()
                    ax3.text(0.35, 0.5, 'Unimodal', fontsize=titlefont)
                else:
                    categories = categories  # unimodal with categories

                    low_exp_group = sum([x == 1 for x in categories])
                    high_exp_group = sum([x == 2 for x in categories])

                    tol = len(categories)
                    if verbosity == True:
                        print(' No. in Low Expression Group: %s (%s%%)' % (low_exp_group, np.round((low_exp_group/tol)*100, 2)), '\n',
                              'No. in High Expression Group: %s (%s%%)' % (
                                  high_exp_group, np.round((high_exp_group/tol)*100, 2)), '\n',
                              'Number of Total:', tol)

                    y = np.round([low_exp_group, high_exp_group], 2)
                    y_p = np.round(
                        [low_exp_group*100/tol, high_exp_group*100/tol], 2)
                    x = np.arange(2)

                    ax3.bar(x, y, width=0.2, align='center', ec='b')
                    ax3.set_xticks(x)
                    ax3.set_xticklabels(['Low', 'High'], fontsize=axisfont)
                    ax3.set_title('Expression Categories of %s' %
                                  ID, fontsize=titlefont)
                    ax3.set_ylabel('Number of Samples', fontsize=axisfont)

                    for i, v in enumerate(y_p):
                        ax3.text(x[i]-0.05, y[i] + max(y)/100,
                                 str(v)+'%', fontsize=axisfont)

            elif isinstance(cutoff, float) | isinstance(cutoff, np.integer):

                if calc_back == False:
                    # determine patients
                    low_exp_group = input_data.columns[(
                        input_data < cutoff).loc[ID]]
                    high_exp_group = input_data.columns[(
                        input_data >= cutoff).loc[ID]]
                else:
                    # determine datapoints
                    low_exp_group = data[data <= cutoff]
                    high_exp_group = data[data > cutoff]

                tol = len(low_exp_group) + len(high_exp_group)
                if verbosity == True:
                    print(' No. in Low Expression Group: %s (%s%%)' % (len(low_exp_group), np.round((len(low_exp_group)/tol)*100, 2)), '\n',
                          'No. in High Expression Group: %s (%s%%)' % (
                              len(high_exp_group), np.round((len(high_exp_group)/tol)*100, 2)), '\n',
                          'Number of Total:', tol)

                y = np.round([len(low_exp_group), len(high_exp_group)], 2)
                y_p = np.round([len(low_exp_group)*100/tol,
                                len(high_exp_group)*100/tol], 2)
                x = np.arange(2)

                ax3.bar(x, y, width=0.2, align='center', ec='b')
                ax3.set_xticks(x)
                ax3.set_xticklabels(['Low', 'High'], fontsize=axisfont)
                ax3.set_title('Expression Categories of %s' %
                              ID, fontsize=titlefont)
                ax3.set_ylabel('Number of Samples', fontsize=axisfont)

                if calc_back == False:
                    # determine datapoints
                    x_lowexp = pd.DataFrame(input_data.loc[:, low_exp_group])
                    x_highexp = pd.DataFrame(input_data.loc[:, high_exp_group])

                    if filt != None:
                        true_posh = x_highexp.columns[(x_highexp > twostd).loc[ID]]

                        categories = [
                            2 if x in true_posh else 1 for x in input_datanormcat]

                        y2 = np.round([len([]), len(true_posh)], 2)

                        ax3.bar(x, y2, width=0.2, align='center',
                                ec='b', facecolor='g')

                        background_d = mpatches.Patch(facecolor='green', ec='black',
                                                      label='True Positive(s): %s samples' % (len(true_posh)))
                        ax3.legend(handles=[background_d])
                    else:
                        categories = [
                            2 if x in x_highexp.columns else 1 for x in input_datanormcat]

                for i, v in enumerate(y_p):
                    ax3.text(x[i]-0.05, y[i] + max(y)/100,
                             str(v)+'%', fontsize=axisfont)

            elif len(cutoff) == 2:
                if calc_back == False:
                    # determine patients
                    low_exp_group = input_data.columns[(
                        input_data < cutoff[0]).loc[ID]]
                    med_exp_group = input_data.columns[(
                        cutoff[0] <= input_data.loc[ID]) & (input_data.loc[ID] < cutoff[1])]
                    high_exp_group = input_data.columns[(
                        input_data >= cutoff[1]).loc[ID]]
                else:
                    # determine datapoints
                    low_exp_group = data[data < cutoff[0]]
                    med_exp_group = data[(cutoff[0] <= data) & (data < cutoff[1])]
                    high_exp_group = data[data >= cutoff[1]]

                tol = len(low_exp_group) + len(high_exp_group)+len(med_exp_group)
                if verbosity == True:
                    print(' No. in Low Expression Group: %s (%s%%)' % (len(low_exp_group), np.round((len(low_exp_group)/tol)*100, 2)), '\n',
                          'No. in Med Expression Group: %s (%s%%)' % (
                              len(med_exp_group), np.round((len(med_exp_group)/tol)*100, 2)), '\n',
                          'No. in High Expression Group: %s (%s%%)' % (
                              len(high_exp_group), np.round((len(high_exp_group)/tol)*100, 2)), '\n',
                          'Number of Total:', tol)

                y = np.round([len(low_exp_group), len(
                    med_exp_group), len(high_exp_group)], 2)
                y_p = np.round([len(low_exp_group)*100/tol, len(med_exp_group)
                                * 100/tol, len(high_exp_group)*100/tol], 2)
                x = np.arange(3)

                ax3.bar(x, y, width=0.2, align='center', ec='b')
                ax3.set_xticks(x)
                ax3.set_xticklabels(['Low', 'Med', 'High'], fontsize=axisfont)
                ax3.set_title('Expression Categories of %s' %
                              ID, fontsize=titlefont)
                ax3.set_ylabel('Number of Samples', fontsize=axisfont)

                if calc_back == False:

                    # determine datapoints
                    x_lowexp = pd.DataFrame(input_data.loc[:, low_exp_group])
                    x_medexp = pd.DataFrame(input_data.loc[:, med_exp_group])
                    x_highexp = pd.DataFrame(input_data.loc[:, high_exp_group])

                    if filt != None:
                        true_posh = x_highexp.columns[(
                            (x_highexp > twostd).loc[ID])]

                        true_posm = x_medexp.columns[((x_medexp > twostd).loc[ID])]
                        categories = [2 if (x in true_posh) | (
                            x in true_posm) else 1 for x in input_datanormcat]

                        y2 = np.round([len([]), len(true_posm), len(true_posh)], 2)

                        ax3.bar(x, y2, width=0.2, align='center',
                                ec='b', facecolor='g')

                        background_d = mpatches.Patch(facecolor='green', ec='black',
                                                      label='True Positive(s) = %s samples' % (len(true_posh)+len(true_posm)))
                        ax3.legend(handles=[background_d])

                    else:
                        categories = [2 if (x in x_highexp.columns) | (
                            x in x_medexp.columns) else 1 for x in input_datanormcat]

                for i, v in enumerate(y_p):
                    ax3.text(x[i]-0.08, y[i]+max(y)/100,
                             str(v)+'%', fontsize=axisfont)

        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------

            cclemon = np.array(cclemon)

            # plot chisquare value
            ax4.tick_params(axis='both', labelsize=axisfont)
            from mpl_toolkits.axes_grid1 import make_axes_locatable

            if condition != 0 and condition != 5:
                # chisquare value graph
                xp = np.arange(0, len(cclemon))
                length = len(cclemon[0])
                if length == 2:

                    color = ['r', 'b']

                    l1 = mlines.Line2D([], [], color=color[0],
                                       label='Fitting w/ 1 n_comp')
                    l2 = mlines.Line2D([], [], color=color[1],
                                       label='Fitting w/ 2 n_comp')
                    arrows = mpatches.FancyArrowPatch(
                        0.1, 0.2, color='r', label='Lower X\N{SUPERSCRIPT TWO} Found')

                    # if the max is not at least twice the initial then dont break axis
                    if np.max(cclemon[np.isfinite(cclemon)])/chis[0] > 2.5:
                        divider = make_axes_locatable(ax4)
                        # break the axis to show lower chi-square
                        ax44 = divider.new_vertical(size="100%", pad=0.1)
                        f.add_axes(ax44)

                        ax4.set_ylim([0, chis[0]*1.1])
                        ax4.spines['top'].set_visible(False)
                        ax44.set_ylim([np.mean(cclemon[np.isfinite(cclemon)]), max(
                            cclemon[np.isfinite(cclemon)])])
                        ax44.tick_params(bottom="off", labelbottom='off')
                        ax44.spines['bottom'].set_visible(False)
                        ax44.legend(handles=[l1, l2, arrows])
                    else:
                        ax4.legend(handles=[l1, l2, arrows])

                    for zz in range(0, length):
                        lemonbags = []
                        for lemons in cclemon:
                            lemonbags.append(lemons[zz])
                        ax4.plot(xp, lemonbags, '-', c=color[zz])
                        # if the max is not at least twice the initial then dont break axis
                        if np.max(cclemon[np.isfinite(cclemon)])/chis[0] > 2.5:
                            ax44.plot(xp, lemonbags, '-', c=color[zz])

                    highlight = [[x, y] for x, y in zip(rem[2:], chis[1:])]

                    for h in highlight:
                        ax4.annotate('(%s)' % h[0], ha='center', fontsize=axisfont, xy=(h[0], h[1]), xytext=(0, 25), textcoords='offset points',
                                     arrowprops=dict(ec='r', lw=4))

                else:
                    highlight = [[x, y] for x, y in zip(rem[2:], chis[1:])]
                    color = ['r', 'b', 'g', 'orange']
                    l1 = mlines.Line2D([], [], color=color[0],
                                       label='Left Tail w/ 1 n_comp')
                    l2 = mlines.Line2D([], [], color=color[1],
                                       label='Left Tail w/ 2 n_comp')
                    r1 = mlines.Line2D([], [], color=color[2],
                                       label='Right Tail w/ 1 n_comp')
                    r2 = mlines.Line2D([], [], color=color[3],
                                       label='Right Tail w/ 2 n_comp')
                    arrows = mpatches.FancyArrowPatch(
                        0.1, 0.2, color='r', label='Lower X\N{SUPERSCRIPT TWO} Found')

                    # if the max is not at least twice the initial then dont break axis
                    if np.max(cclemon[np.isfinite(cclemon)])/chis[0] > 2.5:
                        divider = make_axes_locatable(ax4)
                        # break the axis to show lower chi-square
                        ax44 = divider.new_vertical(size="100%", pad=0.1)
                        f.add_axes(ax44)

                        ax4.set_ylim([0, chis[0]*1.1])
                        ax4.spines['top'].set_visible(False)
                        ax44.set_ylim([np.mean(cclemon[np.isfinite(cclemon)]), max(
                            cclemon[np.isfinite(cclemon)])])
                        ax44.tick_params(bottom="off", labelbottom='off')
                        ax44.spines['bottom'].set_visible(False)
                        ax44.legend(handles=[l1, l2, r1, r2, arrows])
                    else:
                        ax4.legend(handles=[l1, l2, r1, r2, arrows])

                    for zz in range(0, length):
                        lemonbags = []
                        for lemons in cclemon:
                            lemonbags.append(lemons[zz])
                        ax4.plot(xp, lemonbags, '-', c=color[zz])
                        # if the max is not at least twice the initial then dont break axis
                        if np.max(cclemon[np.isfinite(cclemon)])/chis[0] > 2.5:
                            ax44.plot(xp, lemonbags, '-', c=color[zz])

                    for h in highlight:
                        ax4.annotate('(%s)' % h[0], ha='center', fontsize=axisfont, xy=(h[0], h[1]), xytext=(0, 25), textcoords='offset points',
                                     arrowprops=dict(ec='r', lw=4))

                # if the max is not at least twice the initial then dont break axis
                if np.max(cclemon[np.isfinite(cclemon)])/chis[0] > 2.5:
                    # From https://matplotlib.org/examples/pylab_examples/broken_axis.html
                    d = .015  # how big to make the diagonal lines in axes coordinates
                    # arguments to pass to plot, just so we don't keep repeating them
                    kwargs = dict(transform=ax44.transAxes,
                                  color='k', clip_on=False)
                    # top-left diagonal
                    ax44.plot((-d, +d), (-d, +d), **kwargs)
                    ax44.plot((1 - d, 1 + d), (-d, +d), **
                              kwargs)  # top-right diagonal

                    # switch to the bottom axes
                    kwargs.update(transform=ax4.transAxes)
                    # bottom-left diagonal
                    ax4.plot((-d, +d), (1 - d, 1 + d), **kwargs)
                    ax4.plot((1 - d, 1 + d), (1 - d, 1 + d), **
                             kwargs)  # bottom-right diagonal

                    ax44.set_title('Iterative Tail Trimming', fontsize=titlefont)
                    ax4.set_xlabel('Datapoints removed', fontsize=axisfont)
                    ax4.set_ylabel('Chi-square Value', fontsize=axisfont)
                    ax4.yaxis.set_label_coords(1.06, 1)
                    ax44.tick_params(axis='both', labelsize=axisfont)
                else:
                    ax4.set_title('Iterative Tail Trimming', fontsize=titlefont)
                    ax4.set_xlabel('Datapoints removed', fontsize=axisfont)
                    ax4.set_ylabel('Chi-square value', fontsize=axisfont)
            else:
                ax4.text(0.35, 0.5, 'Nothing to see here', fontsize=titlefont)

            if graph_show == True:
                plt.tight_layout()
                plt.show()
            else:
                matplotlib.use('Agg')

        elif graphs == False:
            twostd = meanf+2*stdf

            # get cutoff
            if count == 0:
                if n_comp > 1:
                    xf = group_div
                    classif = 'bimodal'
                else:
                    pass
                    xf = None
                    if categories != None:
                        classif = 'categorical unimodal'
                    else:
                        classif = 'unimodal'
            elif condition == 1:  # remove 1 point means removeing 0 index need to zero index
                xf = datas[remove][0]
                classif = 'unimodal with a tail'
            elif condition == 2:
                xf = np.sort(np.concatenate([group_div, datas[remove]]))
                classif = 'bimodal with a tail'
            elif condition == 3:
                xf = datas[-remove][0]
                classif = 'unimodal with a tail'
            elif condition == 4:
                xf = np.sort(np.concatenate([group_div, datas[-remove]]))
                classif = 'bimodal with a tail'
            elif condition == 5:
                xf = group_div
                classif = 'poorly fitted bimodal'

        # get categories
            cutoff = xf
            try:
                if len(cutoff) == 1:
                    cutoff = cutoff[0]
            except:
                pass

            if sum(cutoff == None) == 1:
                if categories == None:
                    categories = np.ones(input_data.shape[1]).tolist()
                else:
                    categories = categories
            elif isinstance(cutoff, float) | isinstance(cutoff, np.integer):

                # determine patients
                low_exp_group = input_data.columns[(input_data < cutoff).loc[ID]]
                high_exp_group = input_data.columns[(input_data >= cutoff).loc[ID]]

                if calc_back == False:

                    # determine datapoints
                    x_lowexp = pd.DataFrame(input_data.loc[:, low_exp_group])
                    x_highexp = pd.DataFrame(input_data.loc[:, high_exp_group])

                    if filt != None:
                        true_posh = x_highexp.columns[(x_highexp > twostd).loc[ID]]
                        categories = [
                            2 if x in true_posh else 1 for x in input_datanormcat]
                    else:
                        categories = [
                            2 if x in x_highexp.columns else 1 for x in input_datanormcat]

            elif len(cutoff) == 2:
                print(cutoff)
                # determine patients
                low_exp_group = input_data.columns[(
                    input_data < cutoff[0]).loc[ID]]
                med_exp_group = input_data.columns[(
                    cutoff[0] <= input_data.loc[ID]) & (input_data.loc[ID] < cutoff[1])]
                high_exp_group = input_data.columns[(
                    input_data >= cutoff[1]).loc[ID]]

                if calc_back == False:

                    # determine datapoints
                    x_lowexp = pd.DataFrame(input_data.loc[:, low_exp_group])
                    x_medexp = pd.DataFrame(input_data.loc[:, med_exp_group])
                    x_highexp = pd.DataFrame(input_data.loc[:, high_exp_group])

                    if filt != None:
                        true_posh = x_highexp.columns[(x_highexp > twostd).loc[ID]]
                        true_posm = x_medexp.columns[((x_medexp > twostd).loc[ID])]
                        categories = [
                            3 if x in true_posh else 2 if x in true_posm else 1 for x in input_datanormcat]
                    else:
                        categories = [
                            3 if x in x_highexp.columns else 2 if x in x_medexp.columns else 1 for x in input_datanormcat]

        if output_bins == True:
            # send out dataframe for each bin so we know what cells is in eaach bin
            bins_cell = output_bins_func(input_datanorm, bins)
        else:
            bins_cell = None

        if calc_back == True:
            return means, np.sqrt(covars), xf

        else:
            if chis == []:
                return [means, covars, xf], classif, categories, bimodalc, bins_cell, f
            else:
                return [means, covars, xf], classif, categories, chis[-1], bins_cell, f


    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    # v
    ##############################################################################################################################
    ##############################################################################################################################


    def output_bins_func(data, bins):
        index = [[bins[x], bins[x+1]]
                 for x in range(0, len(bins)-1)]  # get all pairs of bins edges

        final_df = [data[(data <= x[1]) & (data >= x[0])].dropna(
            axis=1).columns.tolist() for x in index]

        return pd.DataFrame(final_df)

    # THIS IS IT

    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    # v
    ##############################################################################################################################
    ##############################################################################################################################


    def GMM_modelingt(ID, input_datanormal, log2transform=True, dynamic_binning_s=True, tune_factor=0.99, filt=None, verbosity=False, graphs=True,
                      graph_show=True, farinto=0.01, meanf=0, calc_back=False, stdf=0, chisquaremethod=True, unimodal_categories=True,
                      Single_tail_validation=True, cell_lines=[], cell_line_groupname=[], output_bins=False):

        dynamic_binning_s, tune_factor, filt, verbosity, graphs, graph_show, farinto, meanf, calc_back,  stdf, unimodal_categories, chisquaremethod, Single_tail_validation, output_bins = dynamic_binning_s, tune_factor, filt, verbosity, graphs, graph_show, farinto, meanf, calc_back, stdf, unimodal_categories, chisquaremethod, Single_tail_validation, output_bins

        if calc_back == True:
            idf = input_datanormal.values.flatten()
            input_datanorm = idf[~np.isnan(idf)].reshape(-1, 1)

            if log2transform == True:
                # log2 transform
                # pseudocount so that we do not ignore gene with 0 expression
                input_pseudo = input_datanormal.replace(0, np.nan)
                datalog = np.log2(input_pseudo)
                idflog = datalog.values.flatten()
                data = idflog[~np.isnan(idflog)].reshape(-1, 1)
                input_data = data

        else:
            datanorm = np.array(
                input_datanormal.loc[ID].dropna().values).reshape(-1, 1)
            val = input_datanormal.loc[ID].dropna()

            input_datanorm = pd.DataFrame(
                np.array(val), index=val.index, columns=[ID]).T
            input_datanormcat = pd.DataFrame(np.array(val), index=val.index, columns=[
                                             ID]).T  # used for categories cuz we need NA to count
            if log2transform == True:
                # log2 transform
                # pseudocount so that we do not ignore gene with 0 expression
                input_pseudo = input_datanorm.replace(0, np.nan)
                input_data = pd.DataFrame(
                    np.array(np.log2(val.astype(np.float64))), index=val.index, columns=[ID]).T

                input_pseudo = input_pseudo.loc[ID].dropna(
                ).values.reshape(-1, 1).astype(np.float64)

                data = np.log2(input_pseudo)

        # define verbosity for below input

        # run GMM pipeline
        if log2transform == True:
            weights, categories, group_div, datas, condition, chis, rem, BICs, cclemon, biclemon, realbic, n_comp, n_components, BIC, bimodalc, xs, expected, expectedt, expectedboth, nums, count, remove, dof, means, covars = GMM_pipeline(
                data, log2transform, filt, meanf, stdf, verbosity, farinto, dynamic_binning_s, Single_tail_validation, tune_factor, chisquaremethod, unimodal_categories)
        elif log2transform == False:
            weights, categories, group_div, datas, condition, chis, rem, BICs, cclemon, biclemon, realbic, n_comp, n_components, BIC, bimodalc, xs, expected, expectedt, expectedboth, nums, count, remove, dof, means, covars = GMM_pipeline(
                datanorm, log2transform, filt, meanf, stdf, verbosity, farinto, dynamic_binning_s, Single_tail_validation, tune_factor, chisquaremethod, unimodal_categories)

        # run GMM plot fucntion
        if log2transform == True:
            [means, covars, xf], classif, categories, chi, bins_cell, f = GMM_plot(input_data, data, datanorm, input_datanorm, categories, verbosity, graphs, graph_show, n_comp, xs, expected, expectedt, expectedboth,
                                                                                   ID, log2transform, nums, meanf, stdf, group_div, count, condition,
                                                                                   n_components, BIC, calc_back, input_datanormcat, chis, means,
                                                                                   covars, bimodalc, remove, datas, biclemon, realbic, filt, cell_lines, cclemon,
                                                                                   dof, rem, cell_line_groupname, output_bins)
        elif log2transform == False:
            [means, covars, xf], classif, categories, chi, bins_cell, f = GMM_plot(input_datanorm, datanorm, datanorm, input_datanorm, categories, verbosity, graphs, graph_show, n_comp, xs, expected, expectedt, expectedboth,
                                                                                   ID, log2transform, nums, meanf, stdf, group_div, count, condition,
                                                                                   n_components, BIC, calc_back, input_datanormcat, chis, means,
                                                                                   covars, bimodalc, remove, datas, biclemon, realbic, filt, cell_lines, cclemon,
                                                                                   dof, rem, cell_line_groupname, output_bins)

        return [means, covars, xf, weights], classif, categories, chi, bins_cell, f


    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
    # v
    ##############################################################################################################################
    ##############################################################################################################################

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
