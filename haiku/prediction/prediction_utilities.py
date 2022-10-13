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

from copy import deepcopy
import os
from typing import Any, Dict

import numpy as np
import numpy.typing as npt
from haiku.climate_data.climate import (climMeanNSIDC, loadCESM1Data,
                                        loadNSIDCData)
from haiku.plotting.plotFunctions import (plotAccuracy, plotPredCompare)


def make_predictions(prediction_horizon: int, Vtn: npt.NDArray,
                     Lambda: npt.NDArray) -> npt.NDArray:
    """Make predictions with Koopman model."""
    # forward propagate the eigenvalues (increment power)
    # number of Koopman modes
    ell = Lambda.shape[0]
    # number of snapshots for reconstruction (include initial snapshot)
    m = prediction_horizon

    lambda_pow = np.ones((ell, m), dtype=np.complex_)

    for j in range(1, m):
        lambda_pow[:, j] = Lambda * lambda_pow[:, j-1]

    # create predictions: sum(Koopman mode x coefficient x eigenvalue)
    preds = Vtn @ lambda_pow

    # take only real valued predictions
    preds = preds.real

    # set any negative values to zero
    preds[preds < 0.0] = 0.0

    # set any values greater than 100 to 100
    preds[preds > 100.0] = 100.0

    return preds


def plot_all_comparisons(parameters: Dict[str, Any],
                         full_preds: npt.NDArray,
                         mask: Dict[str, npt.NDArray],
                         lat: npt.NDArray,
                         lon: npt.NDArray,
                         Koopman: Dict[str, npt.NDArray],
                         ensemble_number: str,
                         outputDir: str) -> None:
    """Process data and plot all comparisons.

    Loads data from various files and calls the
    appropriate plotting methods.
    """
    # extract fields from mask dict
    alphaMask = mask["alphaMask"]

    # make a deep copies of the self.parameters dictionary
    # to load NSIDC and CESM1 data
    nsidc_params = deepcopy(parameters)
    cesm1_params = deepcopy(parameters)
    climm_params = deepcopy(parameters)

    # pull NSIDC sea ice concentration data for comparison
    nsidc_params['dataSource'] = 'NSIDC'
    nsidc_params['NSIDCvariableName'] = 'merged'
    nsidc_params['start_dateInt'] = parameters['start_dateInt']
    nsidc_params['stop_dateInt'] = np.minimum(
        parameters['stop_dateInt'], 20201201)
    nsidc_climate_data = loadNSIDCData(nsidc_params)
    NSIDCdata = nsidc_climate_data.data
    NSIDCdateInt = nsidc_climate_data.date_int
    NSIDClat = nsidc_climate_data.latitudes
    NSIDCdata[
        NSIDCdata <= nsidc_params['concentrationToConsiderZero']] = 0
    NSIDCdata[NSIDCdata > 100.0] = 0.0

    # pull CESM1 sea ice concentration data for comparison
    cesm1_params['dataSource'] = 'CESM1'
    cesm1_params['CESM1variableName'] = 'ICEFRAC'
    cesm1_climate_data = loadCESM1Data(ensemble_number, cesm1_params)
    CESM1data = cesm1_climate_data.data
    CESM1dateInt = cesm1_climate_data.date_int
    CESM1data[
        CESM1data <= cesm1_params['concentrationToConsiderZero']] = 0
    CESM1data[CESM1data > 100.0] = 0.0

    # compute the climatological mean to compute prediction skill
    climm_params['dataSource'] = 'NSIDC'
    climm_params['NSIDCvariableName'] = 'merged'
    climMean = climMeanNSIDC(climm_params)
    climMean[
        climMean <= climm_params['concentrationToConsiderZero']] = 0
    climMean[climMean > 100.0] = 0.0

    # trim data based on spatial region
    if parameters['hemisphere'] == 'N':
        inds = np.argwhere(
            NSIDClat >= parameters['north_lat_bound'])
        # flatten list of lists
        inds_list = [ind[0] for ind in inds]
        NSIDCdata = NSIDCdata[:, inds_list, :]
        CESM1data = CESM1data[:, inds_list, :]
        climMean = climMean[:, inds_list, :]
    elif parameters['hemisphere'] == 'S':
        inds = np.argwhere(
            NSIDClat <= parameters['south_lat_bound'])
        # flatten list of lists
        inds_list = [ind[0] for ind in inds]
        NSIDCdata = NSIDCdata[:, inds_list, :]
        CESM1data = CESM1data[:, inds_list, :]
        climMean = climMean[:, inds_list, :]

    NSIDCdata = np.transpose(NSIDCdata)
    CESM1data = np.transpose(CESM1data)
    climMean = np.transpose(climMean)

    # need the climMean to match the interval under consideration
    monthInt_nsidc = (NSIDCdateInt // 100 % 100) - 1
    monthInt_cesm1 = (CESM1dateInt // 100 % 100) - 1
    yearmonth_cesm1 = CESM1dateInt // 100

    climMeanData = climMean[:, :, monthInt_cesm1]

    # apply mask on the data before plotting/computing accuracies
    # full prediction window
    n = full_preds.shape[2]
    for iday in range(n):
        full_preds[:, :, iday] = full_preds[:, :, iday] * alphaMask
        CESM1data[:, :, iday] = CESM1data[:, :, iday] * alphaMask
        climMeanData[:, :, iday] = climMeanData[:, :, iday] * alphaMask
    # window up until 20201201 for NSIDC
    n_truth = NSIDCdata.shape[2]
    for iday in range(n_truth):
        NSIDCdata[:, :, iday] = NSIDCdata[:, :, iday] * alphaMask

    # set nothern latitudes to 100% sea ice concentration
    # to deal with the pole hole
    # this can be removed after the data
    # is interpolated over the pole hole
    top_ind = -7
    full_preds[:, top_ind:, :] = 100.0
    CESM1data[:, top_ind:, :] = 100.0
    NSIDCdata[:, top_ind:, :] = 100.0
    climMeanData[:, top_ind:, :] = 100.0

    # call plotting functions
    nsidc_pred_compare_output_directory = os.path.join(
        outputDir, "nsidc_prediction_compare")
    os.makedirs(nsidc_pred_compare_output_directory, exist_ok=True)
    plotPredCompare(
        parameters,
        Koopman,
        mask,
        lat,
        lon,
        NSIDCdata,
        full_preds,
        nsidc_pred_compare_output_directory)

    cesm1_pred_compare_output_directory = os.path.join(
        outputDir, "cesm1_prediction_compare")
    os.makedirs(cesm1_pred_compare_output_directory, exist_ok=True)
    plotPredCompare(
        parameters,
        Koopman,
        mask,
        lat,
        lon,
        CESM1data,
        full_preds,
        cesm1_pred_compare_output_directory)

    # plot accuracies for sea ice in the month of interest only
    # 8th index corresponds to 9th month (September) for NORTH
    # 2th index corresponds to 3rd month (March) for SOUTH
    plot_month = 8
    month_ind_nsidc = np.squeeze(
        np.argwhere(monthInt_nsidc == plot_month))
    month_ind_cesm1 = np.squeeze(
        np.argwhere(monthInt_cesm1 == plot_month))
    full_predsMonth = full_preds[:, :, month_ind_cesm1]
    NSIDCdataMonth = NSIDCdata[:, :, month_ind_nsidc]
    CESM1dataMonth = CESM1data[:, :, month_ind_cesm1]
    climMeanMonth_ind = plot_month * np.ones(
        (month_ind_cesm1.shape[0],), dtype=int)
    climMeanDataMonth = climMean[:, :, climMeanMonth_ind]
    yearmonth_cesm1 = yearmonth_cesm1[month_ind_cesm1]

    accuracy_output_directory = os.path.join(
        outputDir, "accuracy")
    os.makedirs(accuracy_output_directory, exist_ok=True)
    plotAccuracy(
        parameters,
        Koopman,
        mask,
        lat,
        lon,
        NSIDCdataMonth,
        full_predsMonth,
        CESM1dataMonth,
        climMeanDataMonth,
        yearmonth_cesm1,
        accuracy_output_directory)
