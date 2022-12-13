"""
@classification
UNCLASSIFIED

@itar
ITAR CONTROLLED

@copyright
Copyright BAE Systems
Copyright (c) 2022 AIMdyn Inc.

Please reference the LICENSE.txt file included with this software
package for all terms and conditions related to this software
including license, rights and distribution.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Any, Dict, List
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta
from haiku.climate_data.climate_data import ClimateData, calcRMSE, calcSpatialCorr, calcCumulative


def plotMeanAndAnnualModes(parameters, Koopman, mask, meanModeNums, annualModeNums, lats, lons, outputDir, descStrs):
    """
    Plot mean and annual variation modes.
    """
    # pull variables from Koopman model dictionary
    reomega = Koopman['reomega']
    imomega = Koopman['imomega']
    modes = Koopman['modes']
    mode_imp = Koopman['mode_imp']

    # pull variables from mask dictionary
    alphaMask = mask['alphaMask']

    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']

    numVars = np.shape(modes)[3]

    # plot the selected mean modes
    for imode in range(len(meanModeNums)):
        modeNum = meanModeNums[imode]
        # create title strings with corresponding eigenvalue information
        eig_str = (r'$norm={:.2f}$' '\n'
                   r'$Re(omega)={:.3f} \frac{{1}}{{years}}, Im(omega)={:.3f} \frac{{1}}{{years}},$' '\n'
                   r'$tau_{{decay}}={:.3f} years, freq_{{osc}}={:.3f} \frac{{1}}{{years}}$'
                  ).format( mode_imp[modeNum-1],
                            12*reomega[modeNum-1],
                            12*imomega[modeNum-1],
                            np.absolute(1/(12*reomega[modeNum-1])),
                            np.absolute(12*imomega[modeNum-1]) )

        for ivar in range(numVars):
            title_str = 'Mode {} ({})\n{}'.format(modeNum, descStrs[ivar], eig_str)

            colorbar_label = parameters['colorbar_label'][ivar] #'concentration (percent)'
            clims = parameters['clims'][ivar] #[0, 100]

            # pull variables from mask dictionary
            alphaMask = mask['alphaMask'][ivar]

            # reshape each mode to original snapshot size
            modes_reshaped = np.reshape(modes[:,:,:,ivar].flatten(), (len(lons[ivar]),len(lats[ivar]),-1) , order='F')

            # visualization of a Koopman mode
            fig = plt.figure()
            fig = plotSnapshot(fig, modes_reshaped[:,:,modeNum-1], title_str, colorbar_label, clims, lats[ivar], lons[ivar], alphaMask)
            if(flag_save_plots):
                savename_str = os.path.join(outputDir, 'meanMode_mode{}_{}.png'.format( modeNum, descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
                plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
            if(flag_close_all_plots):
                plt.close(fig)

    # plot the selected annual modes
    for imode in range(len(annualModeNums)):
        modeNum = annualModeNums[imode]
        # create title strings with corresponding eigenvalue information
        eig_str = (r'$norm={:.2f}$' '\n'
                   r'$Re(omega)={:.3f} \frac{{1}}{{years}}, Im(omega)={:.3f} \frac{{1}}{{years}},$' '\n'
                   r'$tau_{{decay}}={:.3f} years, freq_{{osc}}={:.3f} \frac{{1}}{{years}}$'
                  ).format( mode_imp[modeNum-1],
                            12*reomega[modeNum-1],
                            12*imomega[modeNum-1],
                            np.absolute(1/(12*reomega[modeNum-1])),
                            np.absolute(12*imomega[modeNum-1]) )

        for ivar in range(numVars):
            title_str = 'Mode {} ({})\n{}'.format(modeNum, descStrs[ivar], eig_str)

            colorbar_label = parameters['colorbar_label'][ivar] #'concentration (percent)'
            clims = parameters['clims'][ivar] #[0, 100]

            # pull variables from mask dictionary
            alphaMask = mask['alphaMask'][ivar]

            # reshape each mode to original snapshot size
            modes_reshaped = np.reshape(modes[:,:,:,ivar].flatten(), (len(lons[ivar]),len(lats[ivar]),-1) , order='F')

            # visualization of a Koopman mode
            fig = plt.figure()
            fig = plotSnapshot(fig, modes_reshaped[:,:,modeNum-1], title_str, colorbar_label, clims, lats[ivar], lons[ivar], alphaMask)
            if(flag_save_plots):
                savename_str = os.path.join(outputDir, 'annualMode_mode{}_{}.png'.format( modeNum, descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
                plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
            if(flag_close_all_plots):
                plt.close(fig)

    # plot composite mean mode formed by summing the selected modes
    if(len(meanModeNums) >= 1):
        # only take modes with positive imaginary eigenvalues
        inds = np.argwhere( imomega[meanModeNums-1] >= 0 )
        meanModeNums_sel = meanModeNums[inds]
        # summedMeanMode = np.sum(modes_reshaped[:,:,meanModeNums_sel-1], axis=2)
        summedMeanMode = np.sum(modes[:,:,meanModeNums_sel-1,:], axis=2)

        for ivar in range(numVars):
            title_str = 'Sum of mean modes ({})'.format(descStrs[ivar]) # create title string

            colorbar_label = parameters['colorbar_label'][ivar] #'concentration (percent)'
            clims = parameters['clims'][ivar] #[0, 100]

            # pull variables from mask dictionary
            alphaMask = mask['alphaMask'][ivar]

            # reshape each mode to original snapshot size
            summedMeanMode_reshaped = np.reshape(summedMeanMode[:,:,ivar].flatten(), (len(lons[ivar]),len(lats[ivar]),-1) , order='F')

        # visualization of a Koopman mode
        fig = plt.figure()
        fig = plotSnapshot(fig, summedMeanMode_reshaped[:,:,0], title_str, colorbar_label, clims, lats[ivar], lons[ivar], alphaMask)
        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'meanMode_summedModes_{}.png'.format(descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close(fig)

        # plot without specified clim
        fig = plt.figure()
        fig = plotSnapshot(fig, summedMeanMode_reshaped[:,:,0], title_str, colorbar_label, [], lats[ivar], lons[ivar], alphaMask)
        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'meanMode_summedModes_{}_relScale.png'.format(descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close(fig)

    else:
        summedMeanMode = None
        print('No summed mean mode for {}'.format(descStrs))

    # plot composite annual mode formed by summing the selected modes
    if(len(annualModeNums) >= 1):
        # only take modes with positive imaginary eigenvalues
        inds = np.argwhere( imomega[annualModeNums-1] >= 0 )
        annualModeNums_sel = annualModeNums[inds]
        # summedAnnualMode = np.sum(modes_reshaped[:,:,annualModeNums_sel-1], axis=2)
        summedAnnualMode = np.sum(modes[:,:,annualModeNums_sel-1,:], axis=2)

        for ivar in range(numVars):
            title_str = 'Sum of annual modes ({})'.format(descStrs[ivar]) # create title string

            colorbar_label = parameters['colorbar_label'][ivar] #'concentration (percent)'
            clims = parameters['clims'][ivar] #[0, 100]

            # pull variables from mask dictionary
            alphaMask = mask['alphaMask'][ivar]

            # reshape each mode to original snapshot size
            summedAnnualMode_reshaped = np.reshape(summedAnnualMode[:,:,ivar].flatten(), (len(lons[ivar]),len(lats[ivar]),-1) , order='F')

        # visualization of a Koopman mode
        fig = plt.figure()
        fig = plotSnapshot(fig, summedAnnualMode_reshaped[:,:,0], title_str, colorbar_label, clims, lats[ivar], lons[ivar], alphaMask)
        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'annualMode_summedModes_{}.png'.format(descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close(fig)

        # plot without specified clim
        fig = plt.figure()
        fig = plotSnapshot(fig, summedAnnualMode_reshaped[:,:,0], title_str, colorbar_label, [], lats[ivar], lons[ivar], alphaMask)
        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'annualMode_summedModes_{}_relScale.png'.format(descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close(fig)
    else:
        summedAnnualMode = None
        print('No summed annual mode for {}'.format(descStrs))

    return summedMeanMode, summedAnnualMode


def plotEigenvalues(parameters, Koopman, outputDir, descStr):
    """
    Visualizes the Koopman spectrum on the complex plane.
    """
    # pull variables from Koopman model dictionary
    relambda = Koopman['relambda']
    imlambda = Koopman['imlambda']
    reomega = Koopman['reomega']
    imomega = Koopman['imomega']
    mode_imp = Koopman['mode_imp']

    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']
    Nvalues_to_plot = parameters['Nvalues_to_plot']
    flag_eval_scale_plot = parameters['flag_eval_scale_plot']

    # plot the largest N modes
    if(Nvalues_to_plot > len(reomega)):
        Nvalues_to_plot = len(reomega)
    elif(Nvalues_to_plot < 1):
        Nvalues_to_plot = len(reomega)

    # only plot Nvalues_to_plot top eigenvalues
    relambda_t   = relambda[:Nvalues_to_plot]
    imlambda_t   = imlambda[:Nvalues_to_plot]
    reomega_t    = reomega[:Nvalues_to_plot]
    imomega_t    = imomega[:Nvalues_to_plot]
    mode_imp_t = mode_imp[:Nvalues_to_plot]

    # scale scatter plot size/color based on flag_eval_scale_plot
    if(flag_eval_scale_plot):
        plot_sizes  = 1000*np.log( (mode_imp_t-np.min(mode_imp_t) + 1.0)/(np.max(mode_imp_t)-np.min(mode_imp_t)) + 1)
        plot_colors = np.log(mode_imp_t + 1e-14)
    else:
        plot_sizes  = 60
        plot_colors = 'r'

    # create plots using only Nvalues_to_plot top eigenvalues

    # lambda values
    fig = plt.figure()
    unit_disc = plt.Circle((0, 0), 12, fill=False)
    ax = plt.gca()
    ax.add_artist(unit_disc)
    plt.scatter(12*relambda_t, 12*imlambda_t, s=plot_sizes, c=plot_colors )
    plt.xlabel(r'$Re(\lambda)$')
    plt.ylabel(r'$Im(\lambda)$')
    plt.title('Eigenvalues\n{}'.format(descStr))
    plt.grid()
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-13, 13])
    ax.set_ylim([-13, 13])
    if(flag_eval_scale_plot):
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Log Scaled Mode Norm', rotation=90)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'eigenvaluesDT_{}_.png'.format( descStr.replace(' ', '').replace('(','').replace(')','').replace(':','_').replace(',','_') ) )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)

    # omega values
    fig = plt.figure()
    plt.axvline(x=0, c='black', zorder=1)
    plt.axhline(y=0, c='black', zorder=1)
    plt.scatter(12*reomega_t, 12*imomega_t, s=plot_sizes, c=plot_colors, zorder=2 )
    plt.xlabel(r'$1/\tau_{decay} (1/years)$')
    plt.ylabel(r'$freq^{osc} (1/years)$')
    plt.title('Eigenvalues\n{}'.format(descStr))
    plt.grid()
    ax = plt.gca()
    # ax.set_aspect('equal')
    ax.set_xlim([-5, 5])
    ax.set_ylim([-6.3, 6.3])
    if(flag_eval_scale_plot):
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Log Scaled Mode Norm', rotation=90)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'eigenvaluesCT_{}.png'.format( descStr.replace(' ', '').replace('(','').replace(')','').replace(':','_').replace(',','_') ) )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)


def plotEigenmodes(parameters, Koopman, mask, lats, lons, outputDir, descStrs):
    """
    Visualizes the Koopman modes.
    """
    # pull variables from Koopman model dictionary
    reomega = Koopman['reomega']
    imomega = Koopman['imomega']
    modes = Koopman['modes']
    mode_imp = Koopman['mode_imp']

    # pull variables from mask dictionary
    alphaMask = mask['alphaMask']

    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']
    Nmodes_to_plot = parameters['Nmodes_to_plot']

    # plot the largest N modes
    if(Nmodes_to_plot > len(reomega)):
        Nmodes_to_plot = len(reomega)
    elif(Nmodes_to_plot < 1):
        Nmodes_to_plot = len(reomega)

    numVars = parameters['numVars']

    # keep current order of modes
    reomega_s = reomega
    imomega_s = imomega
    modes_s = modes
    mode_imp_s = mode_imp

    # plot the N largest norm modes
    for imode in range(Nmodes_to_plot):
        # create title strings with corresponding eigenvalue information
        eig_str = (r'$norm={:.2f}$' '\n'
                   r'$Re(omega)={:.3f} \frac{{1}}{{years}}, Im(omega)={:.3f} \frac{{1}}{{years}},$' '\n'
                   r'$tau_{{decay}}={:.3f} years, freq_{{osc}}={:.3f} \frac{{1}}{{years}}$'
                  ).format( mode_imp_s[imode],
                            12*reomega_s[imode],
                            12*imomega_s[imode],
                            np.absolute(1/(12*reomega_s[imode])),
                            np.absolute(12*imomega_s[imode]) )

        for ivar in range(numVars):
            # do not plot bulk variables
            if(('bulk' not in descStrs[ivar]) and ('forcing' not in descStrs[ivar])):
                title_str = 'Mode {} ({})\n{}'.format(imode+1, descStrs[ivar], eig_str)
                colorbar_label = parameters['colorbar_label'][ivar]
                clims = parameters['clims'][ivar]

                # pull variables from mask dictionary
                alphaMask = mask['alphaMask'][ivar]

                # reshape each mode to original snapshot size
                mode_reshaped = np.reshape(modes_s[ivar][:,:,imode].flatten(), (len(lons[ivar]),len(lats[ivar]),-1) , order='F')
                # mode_reshaped[mode_reshaped>100.0]=100.0

                # visualization of a Koopman mode
                fig = plt.figure()
                fig = plotSnapshot(fig, mode_reshaped, title_str, colorbar_label, clims, lats[ivar], lons[ivar], alphaMask)
                if(flag_save_plots):
                    savename_str = os.path.join(outputDir, 'mode{}_{}.png'.format( imode+1, descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
                    plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
                if(flag_close_all_plots):
                    plt.close(fig)


def plotEigenfunctions(parameters, Koopman, outputDir):
    """
    Plots heatmap of eigenfunctions over time.
    """
    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']
    start_dateInt = parameters['start_dateInt']
    start_datetime = datetime.strptime(str(start_dateInt), '%Y%m%d')

    Kefun = Koopman['Kefun']

    colorbar_label = 'Eigenfunction Magnitude'

    # visualization of eigenfunctions

    fig = plt.figure()
    plt.imshow(np.log10( np.absolute(Kefun) + 1e-14 ))
    plt.xlabel('Time')
    plt.ylabel('Different Eigenfunctions')
    plt.title("Eigenfunctions of Sorted Modes")
    cbar = plt.colorbar(fraction=0.018, pad=0.04)  # Adjust the colorbar height
    cbar.ax.set_ylabel(colorbar_label, rotation=90)

    # number of date ticks to visualize
    num_ticks = 5
    time_ticks = np.linspace(1, Kefun.shape[1]-1, num=num_ticks)
    time_labels = [ (start_datetime + relativedelta(months=int( time_ticks[i] ))).date() for i in range(num_ticks)]
    ax = plt.gca()
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)

    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'efunctions.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)


def plotModeImportance(parameters, Koopman, outputDir):
    """
    Plots mode importance.
    """
    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']

    mode_imp = Koopman['mode_imp']

    # visualization of the importance metric
    fig = plt.figure()
    plt.plot( np.log10(mode_imp + 1e-14) )
    plt.grid()
    plt.xlim([-10, mode_imp.shape[0]+10])
    plt.xlabel('Different Modes')
    plt.ylabel('Importance Metric Magnitude')
    plt.title("Importance Metric Across All Modes")
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'imp_mags.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)


def plotPredictions(parameters, mask, lat, lon, preds, outputDir):
    """
    Creates a series of snapshots visualizing the predictions.
    """

    # pull variables from mask dictionary
    alphaMask = mask['alphaMask'] #[-1]

    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']

    colorbar_label = 'concentration (percent)'
    clims = [0]

    # figure out how many snapshots
    n = preds.shape[2]

    # use prediction horizon start to initialize date
    # train_stop_dateInt = parameters['train_stop_dateInt']
    train_stop_dateInt = parameters['start_dateInt']
    d0 = datetime.strptime(str(train_stop_dateInt), '%Y%m%d')

    # loop through predictions and call plotSnapshot for each prediction
    for i in range(n):
        # set title as current date
        curr_date = (d0 + relativedelta(months=i))
        title_str = 'Prediction for {}'.format(curr_date.date())

        # visualization of a prediction snapshot
        fig = plt.figure()
        fig = plotSnapshot(fig, preds[:,:,i], title_str, colorbar_label, clims, lat, lon, alphaMask)
        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'pred_{}.png'.format( curr_date.strftime("%Y%m%d") ) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close(fig)


def default_color_dicts(color_dict:dict(), data_keys:List[str],target_key="NSIDC")->dict():
    #currently setting colors and styles.  [KMD, Climate Mean, CESM1, NSIDC] is tacitly assumed for color consistenty, NSIDC can go anywhere in list
    #if plotting dicts are empty or improperly defined, update with defaults
    if color_dict == None or len(color_dict) != len(data_keys):
        default_colors = ['b','g','darkorange','k']
        color_dict = dict()
        for i, key in enumerate(data_keys):
            if key == target_key:
                color_dict[key]='r'
            else:
                color_dict[key]=default_colors[i]
    return(color_dict)
                
def default_style_dicts(style_dict:dict(), data_keys:List[str],target_key="NSIDC")->dict():
    if style_dict == None or len(style_dict) != len(data_keys):
        style_dict =dict()
        for key in data_keys:
            if "K" in key: # select koopman for now
                style_dict[key]='--'
            else:
                style_dict[key]='-'
    return(style_dict)

def plotTimeSeries(x_dict:dict(), y_dict:dict(), target_key:str, y_axis_label:str, output_directory:str, color_dict=None, style_dict=None, errors:dict()={},plot_title:str=None):
    """
    Visualize accuracy of predictions over the testing interval.
    data dictionary has keys as series label and array of climate data values ( monthly values)
    target_key is the key for the training dataset (or the dataset to compare to)
    yearmonth is the formatted date that matches the elements in the ordered list
    errors is a dictionary whose keys match any number (including none) of the keys in data.
    errors values are a list of 3 items: ["error label", -data, +data]. where the error label 
    is a string that describes the error itelf (eg. 1-sigma) and -data and +data are the upper 
    and lower bounds with the same shape as the corresponding item in data. 
    """
    #currently setting colors and styles.  [KMD, Climate Mean, CESM1, NSIDC] is tacitly assumed for color consistenty, NSIDC can go anywhere in list
    #if plotting dicts are empty or improperly defined, update with defaults
    style_dict = default_style_dicts(style_dict, y_dict.keys(),target_key)
    color_dict = default_color_dicts(color_dict, y_dict.keys(),target_key)
    os.makedirs(output_directory, exist_ok=True)
    if plot_title==None:
        plot_title=y_axis_label
        # number of date ticks to visualize
    num_ticks = 11
    time_ticks = np.linspace(np.min(x_dict[target_key]),np.max(x_dict[target_key]), num=num_ticks).astype(int)
    #time_labels = [ (start_datetime + relativedelta(years=time_ticks[i])).year for i in range(num_ticks)]

    # visualization of prediction accuracy over time
    fig = plt.figure(figsize=(8,3))
    for key in y_dict.keys():
        plt.plot(x_dict[key],y_dict[key], style_dict[key], color = color_dict[key], label = key )
    for key in errors.keys():
        plt.fill_between(x_dict[key],
                         errors[key][1],
                         errors[key][2],color=color_dict[key],alpha=0.5)

    plt.grid()
    plt.xlim([np.min(x_dict[target_key]), np.max(x_dict[target_key])])
    plt.xlabel('Years')
    plt.ylabel(y_axis_label)
    plt.title(plot_title)
    ax = plt.gca()
    ax.set_xticks(time_ticks)
    #ax.set_xticklabels(time_labels)
    ax.legend([label for label in y_dict.keys()], loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3)
    plt.locator_params(axis="both", tight=True)
    savename_str = os.path.join(output_directory, str(plot_title)+'.png' )
    plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    plt.close('all')
        
        

def plotRobustness(parameters, mask, lat, lon, data:dict(), target_key:str, yearmonth, outputDir,color_dict=None,style_dict=None):
    """
    Visualize accuracy of predictions over the testing interval.
    data dictionary has keys as series label and array of climate data values ( monthly values)
    target_key is the key for the training dataset (or the dataset to compare to)
    yearmonth is the formatted date that matches the elements in the ordered list
    """

    print("Generating Accuracy Visualizations")

    # pull variables from mask dictionary
    alphaMask = mask['alphaMask']
    maskedLon = mask['maskedLon']
    maskedLat = mask['maskedLat']

    dataSource = parameters['dataSource']
    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']

    #grab
    n_dict = dict()
    for key in data.keys():
        n_dict[key]=data[key].shape[2]

    n_min = np.min(list(n_dict.values()))
    n_max = np.max(list(n_dict.values()))

    start_dateInt = parameters['start_dateInt']
    start_datetime = datetime.strptime(str(start_dateInt), '%Y%m%d')

    # ========================================================= #
    #                          RMSE                             #
    # ========================================================= #

    # compute RMSE of predictions vs. climatological mean vs. CESM1 over time (without considering highest latitudes)
    rmses = dict()
    for key in data.keys():
        if key == target_key:
            continue
        rmses[key]=np.zeros(n_min)
    for iday in range(n_min):
        for key in data.keys():
            if key == target_key:
                continue
            rmses[key][iday] = calcRMSE(data[key][:,:,iday], data[target_key][:,:,iday])

    # number of date ticks to visualize
    num_ticks = 11
    time_ticks = np.linspace(0, n_min-1, num=num_ticks).astype(int)
    time_labels = [ (start_datetime + relativedelta(years=time_ticks[i])).year for i in range(num_ticks)]

    # visualization of prediction accuracy over time
    fig = plt.figure(figsize=(8,3))
    for key in rmses.keys():
        plt.plot(rmses[key], label=key)
    plt.grid()
    plt.xlim([0, n_min-1])
    plt.xlabel('Years')
    plt.ylabel('RMSE')
    plt.title('RMSE of September Sea Ice')
    ax = plt.gca()
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)
    ax.legend([label for label in data.keys()], loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3)
    plt.locator_params(axis="both", tight=True)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'rmse.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)

    # ========================================================= #
    #                   Spatial Correlation                     #
    # ========================================================= #

    # compute spatial correlation of predictions vs. climatological mean over time (without considering highest latitudes)
    scorrs = dict()
    for key in data.keys():
        if key == target_key:
            continue
        scorrs[key]=np.zeros(n_min)
    for iday in range(n_min):
        for key in data.keys():
            if key == target_key:
                continue
            scorrs[key][iday] = calcSpatialCorr(data[key][:,:,iday], data[target_key][:,:,iday])

    # number of date ticks to visualize
    num_ticks = 11
    time_ticks = np.linspace(0, n_min-1, num=num_ticks).astype(int)
    time_labels = [ (start_datetime + relativedelta(years=time_ticks[i])).year for i in range(num_ticks)]

    # visualization of prediction accuracy over time
    fig = plt.figure(figsize=(8,3))
    #    plt.plot(scorrs[key])
    plt.grid()
    plt.xlim([0, n_min-1])
    plt.xlabel('Years')
    plt.ylabel('Spatial Correlation')
    plt.title('Spatial Correlation of September Sea Ice')
    ax = plt.gca()
    for key in scorrs.keys():
        plt.plot( scorrs[key], label = key )
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)
    ax.legend([label for label in data.keys()], loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3)
    plt.locator_params(axis="both", tight=True)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'spatial_correlation.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)

    # ========================================================= #
    #                    Sea Ice Area/Extent                    #
    # ========================================================= #

    cumulatives = dict()
    for key in data.keys():
        cumulatives[key]=np.zeros(n_dict[key])
    for key in data.keys():
        for iday in range(n_dict[key]):
            cumulatives[key][iday] = calcCumulative(data[key][:,:,iday], (lon,lat),"POLAR")

    # explicitly append 2021 sea ice extent of 4.72 million sq km
    # true_cum = np.append(true_cum,  4.72)

    # explicitly append 2013 sea ice extent as an estimate for 2022
    # esti_cum = np.append(true_cum, true_cum[-9])

    # number of date ticks to visualize
    num_ticks = 13
    time_ticks = np.linspace(0, n_max-1, num=num_ticks).astype(int)
    time_labels = [ (start_datetime + relativedelta(years=time_ticks[i])).year for i in range(num_ticks)]

    fig = plt.figure(figsize=(8,3))
    for key in cumulatives.keys():
        plt.plot( cumulatives[key], style_dict[key], color=color_dict[key], label = key )
    plt.grid()
    plt.title('September Sea Ice Area')
    # ax1.plot( esti_cum, ':', color='r')
    # ax1.scatter(n_truth+1, esti_cum[-1], s=70, marker='*', color='r', zorder=999)
    plt.xlim([0, n_max-1])
    ax = plt.gca()
    ax.set_ylim([2.5, 8.0])
    ax.set_xlabel('Years')
    ax.set_ylabel('Sea Ice (million sq km)')
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)
    ax.legend([label for label in data.keys()], loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3)
    #lines, labels = ax.get_legend_handles_labels()
    #ax.legend(lines, labels)
    plt.locator_params(axis="both", tight=True)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'cumulative_extent_predictions.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)



def plotPredCompare(parameters, mask, lat, lon, data, preds, outputDir):
    """
    Creates a series of snapshots visualizing the predictions compared to the true data.
    """

    print("Generating Prediction Comparisons")
    os.makedirs(outputDir, exist_ok=True)
    # pull variables from mask dictionary
    alphaMask = mask['alphaMask']

    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']
    colorbar_label = 'concentration (percent)'

    n = data.shape[2] # figure out how many snapshots

    # use prediction horizon start to initialize date
    d0 = datetime.strptime(str(parameters['start_dateInt']), '%Y%m%d') + relativedelta(months=0)

    # loop through predictions and call plotSnapshot for each prediction
    for i in range(n):
        # set title as current date
        curr_date = (d0 + relativedelta(months=i))
        title_str = '{}'.format(curr_date.date())

        # visualization of a prediction snapshot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,3.5))
        fig.subplots_adjust(wspace=0.3)

        c_1 = ax1.imshow(np.rot90(data[:,:,i], axes=(0,1)), alpha = np.rot90(alphaMask, axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='6', interpolation='nearest')
        # c_1 = ax1.imshow(np.rot90(data[:,:,i], axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='6', interpolation='nearest')
        ax1.set_title("Observational Data")

        c_2 = ax2.imshow(np.rot90(preds[:,:,i], axes=(0,1)), alpha = np.rot90(alphaMask, axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='6', interpolation='nearest')
        # c_2 = ax2.imshow(np.rot90(preds[:,:,i], axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='6', interpolation='nearest')
        ax2.set_title("FKPM Predictions")

        plt.suptitle(title_str)
        fig.supxlabel('Longitude [deg]')
        fig.supylabel('Latitude [deg]')

        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'predcomp_%03i.png' % (i+1) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close("all")


def plotSnapshot(fig, data, title_str, colorbar_label, clims, lat, lon, alphaMask):
    """
    Plots a single snapshot of a mode.
    """

    if(len(clims)==2):
        plt.imshow(np.flipud(np.rot90(data, k=1, axes=(0,1))), alpha = np.flipud(np.rot90(alphaMask, k=1, axes=(0,1))), extent=[lon.min(), lon.max(), lat.min(), lat.max()], vmin=clims[0], vmax=clims[1], aspect='auto', interpolation='nearest')
    else:
        plt.imshow(np.flipud(np.rot90(data, k=1, axes=(0,1))), alpha = np.flipud(np.rot90(alphaMask, k=1, axes=(0,1))), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='auto', interpolation='nearest')

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.title(title_str)
    cbar = plt.colorbar(fraction=0.018, pad=0.04)  # Adjust the colorbar height

    return fig

def plotSnapshotsSideBySide(fig, data1, data2, title_str, colorbar_label, clims, x, y, alphaMask):
    """
    Plots a single snapshot of a mode.
    """

    # visualization of a prediction snapshot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,3.5))
    #fig.subplots_adjust(wspace=0.3)

    c_1 = ax1.imshow(np.rot90(data1, axes=(0,1)), alpha = np.rot90(alphaMask, axes=(0,1)), extent=[x.min(), x.max(), y.min(), y.max()], interpolation='nearest')
    ax1.set_title("Observational Data")

    c_2 = ax2.imshow(np.rot90(data2, axes=(0,1)), alpha = np.rot90(alphaMask, axes=(0,1)), extent=[x.min(), x.max(), y.min(), y.max()], interpolation='nearest')
    ax2.set_title("FKPM Predictions")

    plt.suptitle(title_str)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    #fig.supxlabel('Longitude [deg]')
    #fig.supylabel('Latitude [deg]')
    return fig

def plot_robustness_timeseries(comparison_data:dict(), prediction_data:ClimateData,output_directory,analysis_type:str="coverage",target_key="NSIDC"):
    if analysis_type=="coverage":
        ylabel="Sea Ice Extent"
        ypred = prediction_data.return_coverage()
        prediction_low_band, prediction_high_band = \
            prediction_data.return_coverage(use_errors=True)
        yclim = comparison_data["climatological mean"].return_coverage()
        ynsidc = comparison_data["NSIDC"].return_coverage()
    elif analysis_type=="spatial_correlation":
        ylabel = "spatial correlation"
        ypred = prediction_data.return_spatial_correlation(comparison_data[target_key].data)
        prediction_low_band, prediction_high_band = \
            prediction_data.return_spatial_correlation(comparison_data[target_key].data,use_errors=True)
        yclim = comparison_data["climatological mean"].return_spatial_correlation(comparison_data[target_key].data)
        ynsidc = comparison_data["NSIDC"].return_spatial_correlation(comparison_data[target_key].data)
    elif analysis_type=="rmse":
        ylabel="RMSE"
        ypred = prediction_data.return_rmse(comparison_data[target_key].data)
        prediction_low_band, prediction_high_band = \
            prediction_data.return_rmse(comparison_data[target_key].data,use_errors=True)
        yclim = comparison_data["climatological mean"].return_rmse(comparison_data[target_key].data)
        ynsidc = comparison_data["NSIDC"].return_rmse(comparison_data[target_key].data)

    xpred = prediction_data.date_int // 10000
    xnsidc = comparison_data["NSIDC"].date_int // 10000

    for imonth in range(12):
        pred_filter = ((prediction_data.date_int // 100) %100) -1 == imonth
        nsidc_filter =((comparison_data["NSIDC"].date_int // 100) %100) -1 == imonth
        plotTimeSeries({"KMD Prediction":xpred[pred_filter],
                        "Climatological Mean":xnsidc[nsidc_filter],
                        "NSIDC":xnsidc[nsidc_filter]},
                       {"KMD Prediction":ypred[pred_filter],
                        "Climatological Mean":yclim[nsidc_filter],
                        "NSIDC":ynsidc[nsidc_filter]},
                       "NSIDC",
                       ylabel,
                       output_directory,
                       color_dict={},
                       style_dict={},
                       errors={"KMD Prediction":["2-sigma",
                                                 prediction_low_band[pred_filter],
                                                 prediction_high_band[pred_filter]]},
                       plot_title="Month_"+str(imonth+1)+"_"+ylabel.replace(" ","_"))



def plot_timeseries(comparison_data:dict(), prediction_data:ClimateData,output_directory:str,analysis_type:str="coverage",target_key="NSIDC"):

    data_dict={}
    xpred = prediction_data.date_int // 10000
    xnsidc = comparison_data["NSIDC"].date_int // 10000
    yvals = {}
    if analysis_type=="coverage":
        ylabel="Sea Ice Extent"
        yvals["KMD Prediction"] = prediction_data.return_coverage()
        for key in comparison_data.keys():
            yvals[key] = comparison_data[key].return_coverage()
    elif analysis_type=="spatial_correlation":
        ylabel = "spatial correlation"
        yvals["KMD Prediction"] = prediction_data.return_spatial_correlation(comparison_data[target_key].data)
        for key in comparison_data.keys():
            yvals[key] = comparison_data[key].return_spatial_correlation(comparison_data[target_key].data)
    elif analysis_type=="rmse":
        ylabel="RMSE"
        yvals["KMD Prediction"] = prediction_data.return_rmse(comparison_data[target_key].data)
        for key in comparison_data.keys():
            yvals[key]= comparison_data[key].return_rmse(comparison_data[target_key].data)

    for imonth in range(12):
        pred_filter = ((prediction_data.date_int // 100) %100) -1 == imonth
        nsidc_filter =((comparison_data["NSIDC"].date_int // 100) %100) -1 == imonth
        tempx = {"KMD Prediction":xpred[pred_filter]}
        tempy = {"KMD Prediction":yvals["KMD Prediction"][pred_filter]}
        for key in comparison_data.keys():
            tempx[key] = xnsidc[nsidc_filter]
            tempy[key] = yvals[key][pred_filter] 
        plotTimeSeries(tempx,
                       tempy,
                       "NSIDC",
                       ylabel,
                       output_directory,
                       color_dict={},
                       style_dict={},
                       plot_title="Month_"+str(imonth+1)+"_"+ylabel.replace(" ","_"))


    '''
    ===============================================================================
    Copyright (c) 2022 AIMdyn Inc.
    All right reserved.

    3-Clause BSD License

    Additional copyrights may follow.


    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.
    3. Neither the name of the copyright holder nor the names of its contributors
       may be used to endorse or promote products derived from this software
       without specific prior written permission.

    The copyright holder provides no reassurances that the source code
    provided does not infringe any patent, copyright, or any other
    intellectual property rights of third parties.  The copyright holder
    disclaims any liability to any recipient for claims brought against
    recipient by any third party for infringement of that parties
    intellectual property rights.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    ===============================================================================
    '''
