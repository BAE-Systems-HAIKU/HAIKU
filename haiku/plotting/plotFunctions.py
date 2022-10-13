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
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta


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
            modes_reshaped = np.reshape(modes[:,:,:,ivar].flatten(), (len(lons[0]),len(lats[0]),-1) , order='F')

            # visualization of a Koopman mode
            fig = plt.figure()
            fig = plotSnapshot(fig, modes_reshaped[:,:,modeNum-1], title_str, colorbar_label, clims, lats[0], lons[0], alphaMask)
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
            modes_reshaped = np.reshape(modes[:,:,:,ivar].flatten(), (len(lons[0]),len(lats[0]),-1) , order='F')

            # visualization of a Koopman mode
            fig = plt.figure()
            fig = plotSnapshot(fig, modes_reshaped[:,:,modeNum-1], title_str, colorbar_label, clims, lats[0], lons[0], alphaMask)
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
            summedMeanMode_reshaped = np.reshape(summedMeanMode[:,:,ivar].flatten(), (len(lons[0]),len(lats[0]),-1) , order='F')

        # visualization of a Koopman mode
        fig = plt.figure()
        fig = plotSnapshot(fig, summedMeanMode_reshaped[:,:,0], title_str, colorbar_label, clims, lats[0], lons[0], alphaMask)
        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'meanMode_summedModes_{}.png'.format(descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close(fig)

        # plot without specified clim
        fig = plt.figure()
        fig = plotSnapshot(fig, summedMeanMode_reshaped[:,:,0], title_str, colorbar_label, [], lats[0], lons[0], alphaMask)
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
            summedAnnualMode_reshaped = np.reshape(summedAnnualMode[:,:,ivar].flatten(), (len(lons[0]),len(lats[0]),-1) , order='F')

        # visualization of a Koopman mode
        fig = plt.figure()
        fig = plotSnapshot(fig, summedAnnualMode_reshaped[:,:,0], title_str, colorbar_label, clims, lats[0], lons[0], alphaMask)
        if(flag_save_plots):
            savename_str = os.path.join(outputDir, 'annualMode_summedModes_{}.png'.format(descStrs[ivar].replace(' ', '').replace('(','').replace(')','') ) )
            plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
        if(flag_close_all_plots):
            plt.close(fig)

        # plot without specified clim
        fig = plt.figure()
        fig = plotSnapshot(fig, summedAnnualMode_reshaped[:,:,0], title_str, colorbar_label, [], lats[0], lons[0], alphaMask)
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
                mode_reshaped = np.reshape(modes_s[ivar][:,:,imode].flatten(), (len(lons[0]),len(lats[0]),-1) , order='F')
                # mode_reshaped[mode_reshaped>100.0]=100.0

                # visualization of a Koopman mode
                fig = plt.figure()
                fig = plotSnapshot(fig, mode_reshaped, title_str, colorbar_label, clims, lats[0], lons[0], alphaMask)
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


def calcRMSE(pred, true):
    """
    Calculate RMSE given predicted and true values. Both inputs are 2D gridded points.
    """
    rmse = np.sqrt( np.mean( np.square( pred - true ), axis=(0,1) ) )
    return rmse


def calcSpatialCorr(pred, true):
    """
    Calculate spatial correlation given predicted and true values. Both inputs are 2D gridded points.

    This needs some more refining before results are shown.
    """
    scorr = np.corrcoef( pred.flatten(), true.flatten() )[0,1]
    return scorr


def calcCumulative(pred, lat, yearmonth):
    """
    Calculate cumulative given 2D gridded points.
    """

    num_lat = pred.shape[1] #54  # number of latitudes
    num_lon = pred.shape[0] #288 # number of longitudes

    # lat_inc = 0.9375 # 180/192
    lat_inc = 0.94240837696 # 180/(192-1)
    lon_inc = 1.25   # 360/288

    # # latitudes are from -90 to 90 with lat_inc degree as the interval
    # lat = np.zeros(num_lat)
    # for ilat in range(num_lat):
    #     lat[ilat] = -90.0 + lat_inc*ilat
    # # print(lat) # to delete

    earth_rad = 6371.0  # radius of the Earth in km
    earth_deg = (2*np.pi*earth_rad)/360  # 1 degree of the Earth in km
    grid_area = (earth_deg*lat_inc)*(earth_deg*lon_inc)  # area of each lat_inc by lon_inc grid in sq km

    x = 0
    y = 0
    miss_val = -9.99e8  # missing value over land grids
    conc_cutoff = 15    # only consider grids with more that 15% coverage of sea ice
    for ilon in range(num_lon):
        for ilat in range(num_lat):
            if( (pred[ilon,ilat]>=conc_cutoff) and (pred[ilon,ilat]!=miss_val) ):
                x = x + pred[ilon,ilat] * np.cos( (lat[ilat]-lat_inc) /180*np.pi) * grid_area
                # x = x + 100 * np.cos( (lat[ilat]-lat_inc) /180*np.pi) * grid_area
                y = y + 1.0 * np.cos( (lat[ilat]-lat_inc) /180*np.pi) * grid_area

    # 100 is to convert from percent
    # 1000000 is to convert to million sq km.
    x = x / 100.0 / 1000000.0
    y = y / 1000000.0

    return x # sea ice area
    # return y # sea ice extent


def plotAccuracy(parameters, Koopman, mask, lat, lon, NSIDCdata, preds, CESM1data, climMeanData, yearmonth, outputDir):
    """
    Visualize accuracy of predictions over the testing interval.
    """

    print("Generating Accuracy Visualizations")

    # pull variables from mask dictionary
    alphaMask = mask['alphaMask']
    maskedLon = mask['maskedLon']
    maskedLat = mask['maskedLat']

    dataSource = parameters['dataSource']
    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']

    n_pred = preds.shape[2] # number of snapshots
    n_truth = NSIDCdata.shape[2] # number of NSIDC snapshots

    n = np.minimum(n_pred, n_truth)

    pred_horizon = parameters['pred_horizon']
    start_dateInt = parameters['start_dateInt']
    start_datetime = datetime.strptime(str(start_dateInt), '%Y%m%d')

    # ========================================================= #
    #                          RMSE                             #
    # ========================================================= #

    # compute RMSE of predictions vs. climatological mean vs. CESM1 over time (without considering highest latitudes)
    pred_mse = np.zeros(n)
    ref_mse  = np.zeros(n)
    clim_mse = np.zeros(n)
    for iday in range(n):
        pred_mse[iday] = calcRMSE(preds[:,:,iday],        NSIDCdata[:,:,iday])
        ref_mse[iday]  = calcRMSE(CESM1data[:,:,iday],    NSIDCdata[:,:,iday])
        clim_mse[iday] = calcRMSE(climMeanData[:,:,iday], NSIDCdata[:,:,iday])

    # number of date ticks to visualize
    num_ticks = 11
    time_ticks = np.linspace(0, n-1, num=num_ticks).astype(int)
    time_labels = [ (start_datetime + relativedelta(years=time_ticks[i])).year for i in range(num_ticks)]

    # visualization of prediction accuracy over time
    fig = plt.figure(figsize=(8,3))
    plt.plot( pred_mse )
    plt.plot( ref_mse  )
    plt.plot( clim_mse )
    plt.grid()
    plt.xlim([0, n-1])
    plt.xlabel('Years')
    plt.ylabel('RMSE')
    plt.title('RMSE of September Sea Ice')
    ax = plt.gca()
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)
    ax.legend(['KMD Prediction', 'CESM1', 'Climatological Mean'], loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3)
    plt.locator_params(axis="both", tight=True)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'skill_rmse.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)

    # rel_skill = 1 - pred_mse/ref_mse

    # # visualization of prediction skill over time
    # fig = plt.figure(figsize=(8,3))
    # plt.plot(rel_skill)
    # plt.grid()
    # plt.xlim([0, n-1])
    # plt.ylim([-1, 1])
    # plt.xlabel('Months in Future')
    # plt.ylabel('Relative Skill')
    # plt.title('KMD Prediction Skill')
    # ax = plt.gca()
    # # ax.set_xticks(time_ticks)
    # # ax.set_xticklabels(time_labels)
    # ax.legend(['KMD Skill','Prediction Start'], loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.locator_params(axis="both", integer=True, tight=True)
    # if(flag_save_plots):
    #     savename_str = os.path.join(outputDir, 'skill_rel.png' )
    #     plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    # if(flag_close_all_plots):
    #     plt.close(fig)

    # ========================================================= #
    #                   Spatial Correlation                     #
    # ========================================================= #

    # compute spatial correlation of predictions vs. climatological mean over time (without considering highest latitudes)
    pred_scorr = np.zeros(n)
    ref_scorr  = np.zeros(n)
    clim_scorr = np.zeros(n)
    for iday in range(n):
        pred_scorr[iday] = calcSpatialCorr(preds[:,:,iday],        NSIDCdata[:,:,iday])
        ref_scorr[iday]  = calcSpatialCorr(CESM1data[:,:,iday],    NSIDCdata[:,:,iday])
        clim_scorr[iday] = calcSpatialCorr(climMeanData[:,:,iday], NSIDCdata[:,:,iday])

    # number of date ticks to visualize
    num_ticks = 11
    time_ticks = np.linspace(0, n-1, num=num_ticks).astype(int)
    time_labels = [ (start_datetime + relativedelta(years=time_ticks[i])).year for i in range(num_ticks)]

    # visualization of prediction accuracy over time
    fig = plt.figure(figsize=(8,3))
    plt.plot( pred_scorr )
    plt.plot( ref_scorr  )
    plt.plot( clim_scorr )
    plt.grid()
    plt.xlim([0, n-1])
    plt.xlabel('Years')
    plt.ylabel('Spatial Correlation')
    plt.title('Spatial Correlation of September Sea Ice')
    ax = plt.gca()
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)
    ax.legend(['KMD Prediction', 'CESM1', 'Climatological Mean'], loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3)
    plt.locator_params(axis="both", tight=True)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'skill_scorr.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)

    # ========================================================= #
    #                    Sea Ice Area/Extent                    #
    # ========================================================= #

    if(dataSource == 'NSIDC'):
        pred_nsidc_flag = True
    elif(dataSource == 'CESM1'):
        pred_nsidc_flag = False

    # compute cumulative
    pred_cum = np.zeros(n_pred)
    cesm_cum = np.zeros(n_pred)
    true_cum = np.zeros(n_truth)
    clim_cum = np.zeros(n_pred)
    for iday in range(n_pred):
        pred_cum[iday] = calcCumulative(preds[:,:,iday],        lat, yearmonth[iday])
        cesm_cum[iday] = calcCumulative(CESM1data[:,:,iday],    lat, yearmonth[iday])
        clim_cum[iday] = calcCumulative(climMeanData[:,:,iday], lat, yearmonth[iday])
    for iday in range(n_truth):
        true_cum[iday] = calcCumulative(NSIDCdata[:,:,iday],    lat, yearmonth[iday])

    # explicitly append 2021 sea ice extent of 4.72 million sq km
    # true_cum = np.append(true_cum,  4.72)

    # explicitly append 2013 sea ice extent as an estimate for 2022
    # esti_cum = np.append(true_cum, true_cum[-9])

    # number of date ticks to visualize
    num_ticks = 13
    time_ticks = np.linspace(0, n_pred-1, num=num_ticks).astype(int)
    time_labels = [ (start_datetime + relativedelta(years=time_ticks[i])).year for i in range(num_ticks)]

    fig, ax1 = plt.subplots(figsize=(8,3))
    plt.grid()
    ax1.plot( clim_cum, '-', color='g',          label = 'Climatological Mean' )
    ax1.plot( true_cum, '-', color='r',          label = 'NSIDC' )
    # ax1.plot( esti_cum, ':', color='r')
    # ax1.scatter(n_truth+1, esti_cum[-1], s=70, marker='*', color='r', zorder=999)
    ax2 = ax1.twinx()
    ax2.plot( cesm_cum, '-', color='darkorange', label = 'CESM1' )
    if(pred_nsidc_flag):
        ax1.plot( pred_cum, '--', color='b',      label = 'KMD Prediction' )
    else:
        ax2.plot( pred_cum, '--', color='b',      label = 'KMD Prediction' )
    plt.xlim([0, n_pred-1])
    ax1.set_ylim([2.5, 8.0])
    ax2.set_ylim([2.3, 7.8])
    ax1.set_xlabel('Years')
    ax1.set_ylabel('NSIDC Sea Ice (million sq km)', color='r')
    ax2.set_ylabel('CESM1 Sea Ice (million sq km)', color='darkorange')
    plt.title('September Sea Ice Area')
    ax = plt.gca()
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=4)
    plt.locator_params(axis="both", tight=True)
    if(flag_save_plots):
        savename_str = os.path.join(outputDir, 'skill_cum.png' )
        plt.savefig(savename_str, bbox_inches='tight', pad_inches=.1)
    if(flag_close_all_plots):
        plt.close(fig)


def plotPredCompare(parameters, Koopman, mask, lat, lon, NSIDCdata, preds, outputDir):
    """
    Creates a series of snapshots visualizing the predictions compared to the true data.
    """

    print("Generating Prediction Comparisons")

    # pull variables from mask dictionary
    alphaMask = mask['alphaMask']

    flag_save_plots = parameters['flag_save_plots']
    flag_close_all_plots = parameters['flag_close_all_plots']
    colorbar_label = 'concentration (percent)'

    n = NSIDCdata.shape[2] # figure out how many snapshots

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

        c_1 = ax1.imshow(np.rot90(NSIDCdata[:,:,i], axes=(0,1)), alpha = np.rot90(alphaMask, axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='6', interpolation='nearest')
        # c_1 = ax1.imshow(np.rot90(NSIDCdata[:,:,i], axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='6', interpolation='nearest')
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
            plt.close(fig)


def plotSnapshot(fig, data, title_str, colorbar_label, clims, lat, lon, alphaMask):
    """
    Plots a single snapshot of a mode.
    """

    if(len(clims)==2):
        plt.imshow(np.rot90(data, axes=(0,1)), alpha = np.rot90(alphaMask, axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], vmin=clims[0], vmax=clims[1], aspect='auto', interpolation='nearest')
    else:
        plt.imshow(np.rot90(data, axes=(0,1)), alpha = np.rot90(alphaMask, axes=(0,1)), extent=[lon.min(), lon.max(), lat.min(), lat.max()], aspect='auto', interpolation='nearest')

    plt.xlabel('Longitude [deg]')
    plt.ylabel('Latitude [deg]')
    plt.title(title_str)
    cbar = plt.colorbar(fraction=0.018, pad=0.04)  # Adjust the colorbar height
    cbar.ax.set_ylabel(colorbar_label, rotation=90)

    return fig

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