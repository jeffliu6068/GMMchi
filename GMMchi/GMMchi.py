from .GMMchi_func import *
form .dynamic_binning import *
from .solve_gaussians import *

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
