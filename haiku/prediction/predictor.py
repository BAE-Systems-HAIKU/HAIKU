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

import os
from typing import Any, Dict, List

import numpy as np
from haiku.plotting.plotFunctions import (plotEigenvalues,
                                          plotPredictions)
from haiku.prediction.prediction_utilities import (make_predictions,
                                                   plot_all_comparisons)
from haiku.training.dmdR4Haiku import SLS_Spectral_Reconstruct_W_NE
from haiku.training.koopman_model import KoopmanModel


class Predictor:

    def __init__(self, parameters: Dict[str, Any]):
        """Construct a Predictor object."""
        # TODO: Extract parameters into a configuration object
        self.parameters = parameters

    def predict(self, model: KoopmanModel, end_date: int,
                predict_modes: List[int], output_directory: str):
        """Initiate prediction processes with loaded config."""
        self.parameters["stop_dateInt"] = end_date
        self.parameters["predictModes"] = predict_modes

        print('Output directory: {}'.format(output_directory))

        self.parameters['dataSource'] = model.data_source
        self.parameters['NSIDCvariableName'] = model.observable_name
        self.parameters['hemisphere'] = model.hemisphere
        self.parameters['dataFreq'] = model.data_frequency
        self.parameters['runType'] = '40memberLargeEnsemble'
        self.parameters['CESM1variableName'] = 'ICEFRAC'
        self.parameters['ensembleNums'] =\
            list(range(1, 36))+list(range(101, 106))

        # pull training start and end from koopman model
        self.parameters['train_start_dateInt'] =\
            model.training_start_date
        self.parameters['train_stop_dateInt'] =\
            model.training_stop_date

        # select time range for prediction
        self.parameters['start_dateInt'] =\
            self.parameters['train_start_dateInt']

        # determine training horizon
        train_year_diff = int(str(
            self.parameters['train_stop_dateInt'])[0:4]) - \
            int(str(self.parameters['train_start_dateInt'])[0:4])
        train_month_diff = int(str(
            self.parameters['train_stop_dateInt'])[4:6]) - \
            int(str(self.parameters['train_start_dateInt'])[4:6])
        # add 1 here to include the stop date
        self.parameters['train_horizon'] =\
            train_year_diff * 12 + train_month_diff + 1

        # determine prediction horizon
        pred_year_diff = int(str(self.parameters['stop_dateInt'])[0:4]) - \
            int(str(self.parameters['start_dateInt'])[0:4])
        pred_month_diff = int(str(self.parameters['stop_dateInt'])[4:6]) - \
            int(str(self.parameters['start_dateInt'])[4:6])
        # add 1 here to include the stop date
        self.parameters['pred_horizon'] =\
            pred_year_diff * 12 + pred_month_diff + 1

        # ========================================================= #
        # ========================================================= #

        print('\nKoopman Mode Decomposition Based Prediction on Sea Ice Data')
        print('Predictions for {} data using {} modes'.format(
            self.parameters['dataFreq'], len(self.parameters['predictModes'])))
        print('Training Interval:   {} to {} ({} Snapshots)'.format(
            self.parameters['train_start_dateInt'],
            self.parameters['train_stop_dateInt'],
            self.parameters['train_horizon']))
        print('Prediction Interval: {} to {} ({} Snapshots)'.format(
            self.parameters['start_dateInt'],
            self.parameters['stop_dateInt'],
            self.parameters['pred_horizon']))

        # load saved self.parameters
        self.parameters['concentrationToConsiderZero'] =\
            model.concentration_to_consider_zero
        self.parameters['roundConc'] = model.round_concentration
        self.parameters['north_lat_bound'] = model.north_lat_bound
        self.parameters['south_lat_bound'] = model.south_lat_bound
        self.parameters['numSeaIceObsv'] =\
            model.number_sea_ice_observables

        # load all necessary data
        # training data
        data = model.training_data
        # latitude coordinates
        lat = model.latitudes
        # longitude coordinates
        lon = model.longitudes
        self.parameters['nlat'] = len(lat)
        self.parameters['nlon'] = len(lon)

        # load mask
        mask = {}
        mask['maskedLat'] = model.masked_latitudes
        mask['maskedLon'] = model.masked_longitudes

        # load Koopman model
        Koopman = {}
        Koopman['Vtn_real'] = model.vtn_real
        Koopman['Vtn_imag'] = model.vtn_imaginary
        Koopman['relambda'] = model.relambda
        Koopman['imlambda'] = model.imlambda
        Koopman['mode_imp'] = model.mode_imp

        # trim Koopman modes and eigenvalues based on predictModes
        Koopman['Vtn_real'] =\
            Koopman['Vtn_real'][:, self.parameters['predictModes']]
        Koopman['Vtn_imag'] =\
            Koopman['Vtn_imag'][:, self.parameters['predictModes']]
        Koopman['relambda'] =\
            Koopman['relambda'][self.parameters['predictModes']]
        Koopman['imlambda'] =\
            Koopman['imlambda'][self.parameters['predictModes']]
        Koopman['mode_imp'] =\
            Koopman['mode_imp'][self.parameters['predictModes']]

        # complex valued Koopman modes
        Koopman['Vtn'] = Koopman['Vtn_real'] + Koopman['Vtn_imag']*1j
        # complex valued Koopman eigenvalues
        Koopman['Lambda'] = Koopman['relambda'] + Koopman['imlambda']*1j

        m = data.shape[1]

        # set weights for spectral reconstruction
        W = np.ones((m, 1))

        # save coefficients of expansion and
        # scaled Koopman modes for predictions
        Koopman['coefficient'], Koopman['Vtn'] = \
            SLS_Spectral_Reconstruct_W_NE(
                Koopman['Vtn'],
                Koopman['Lambda'],
                data.astype(complex),
                W)

        # generate KMD-based predictions
        preds = make_predictions(
            self.parameters["pred_horizon"],
            Koopman["Vtn"],
            Koopman["Lambda"])
        preds[preds <= self.parameters['concentrationToConsiderZero']] = 0
        preds[preds > 100.0] = 100.0

        np.savetxt(
            os.path.join(output_directory, 'KMD_predictions.csv'),
            preds[:self.parameters['numSeaIceObsv'], :],
            delimiter=",")
        n_pred = preds.shape[1]

        # reshape predictions to original longitude-latitude coordinates
        full_preds = np.zeros(
            (self.parameters['nlon'], self.parameters['nlat'], n_pred))
        for ipred in range(n_pred):
            full_preds[mask['maskedLon'], mask['maskedLat'], ipred] =\
                preds[:self.parameters['numSeaIceObsv'], ipred]

        # apply mask for plotting
        alphaMask = np.zeros(
            (self.parameters['nlon'], self.parameters['nlat']))
        alphaMask[mask['maskedLon'], mask['maskedLat']] = 1.0
        mask['alphaMask'] = alphaMask

        # visualize Koopman eigenvalues used for prediction
        dt = 1.0
        omega = np.log(Koopman['Lambda'])/dt
        Koopman['reomega'] = omega.real*dt
        Koopman['imomega'] = (omega.imag)/(2*np.pi)

        # visualize predictions
        if self.parameters['flag_plot_predictions']:
            prediction_output_directory = os.path.join(
                output_directory, "predictions")
            os.makedirs(prediction_output_directory, exist_ok=True)
            plotPredictions(
                self.parameters,
                mask,
                lat,
                lon,
                full_preds,
                prediction_output_directory)

        # visualize comparisons (need to pull NSIDC data)
        if self.parameters['flag_plot_comparisons']:
            comparisons_output_directory = os.path.join(
                output_directory, "comparisons")
            os.makedirs(comparisons_output_directory, exist_ok=True)

            # NSIDC data doesn't have an ensemble num
            # so we set it to 1 for compatability with
            # CESM1 data (first CESM1 ensemble will be loaded)
            if model.data_source == "CESM1":
                ensemble_number = model.ensemble_number
            else:
                ensemble_number = 1

            plot_all_comparisons(
                self.parameters,
                full_preds,
                mask,
                lat,
                lon,
                Koopman,
                ensemble_number,
                comparisons_output_directory
            )
