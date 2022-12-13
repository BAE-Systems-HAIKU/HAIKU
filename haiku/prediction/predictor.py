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
from haiku.plotting.plotFunctions import plotPredictions
from haiku.climate_data.climate_data import ClimateData
from haiku.prediction.prediction_utilities import make_predictions
from haiku.training.dmdR4Haiku import SLS_Spectral_Reconstruct_W_NE
from haiku.training.koopman_model import KoopmanModel

import datetime
from haiku.utilities.mask import create_mask_dict
from haiku.utilities.date_utilities import (n_dates_in_range,
                                            date_array_for_range)

class Predictor:

    def __init__(self, koopman_model: KoopmanModel):
        """Construct a Predictor object."""
        # TODO: Extract parameters into a configuration object
        self.parameters = dict()
        self.model = koopman_model

        # determine training horizon
        #TODO: shold be defined using self.model.data_frequency;
        #but for now, we only have 1 option, Monthly, so this is okay
        self.parameters['train_horizon'] = n_dates_in_range(self.model.training_start_date,
                                                            self.model.training_stop_date)
        #This seems to be the only configuration parameter not already in the koopman model itself
        #If the predict function is going to live on an object and not as a standalone function,
        #it should be instantiated for a single koopman model at a time
        
    def predict(self, start_date: int, end_date: int, predict_modes: List[int],
                output_directory: str)->ClimateData:
        """Initiate prediction processes with loaded config.
        saves prediction output to .csv in output_directory and returns ClimateData object 
        of the predictions"""
        # select time range for prediction
        if self.model.training_start_date > start_date:
            print("WARNING, start date:", start_date,"is before start date of model training:",self.model.training_start_date)
            print("Overriding start date with training start date:",self.model.training_start_date,"for prediction")
            start_date = self.model.training_start_date

        os.makedirs(output_directory, exist_ok=True)
        print('Output directory: {}'.format(output_directory))


        # determine prediction horizon
        pred_horizon = n_dates_in_range(start_date,end_date)
        
        # ========================================================= #
        # ========================================================= #

        print('\nKoopman Mode Decomposition Based Prediction on Sea Ice Data')
        print('Predictions for {} data using {} modes'.format(
            self.model.data_frequency, len(predict_modes)))
        print('Training Interval:   {} to {} ({} Snapshots)'.format(
            self.model.training_start_date,
            self.model.training_stop_date,
            self.parameters['train_horizon']))
        print('Prediction Interval: {} to {} ({} Snapshots)'.format(
            start_date,
            end_date,
            pred_horizon))

        # load all necessary data
        # training data
                        
        self.parameters['nlat'] = len(self.model.latitudes)
        self.parameters['nlon'] = len(self.model.longitudes)

        # load mask
        #TODO: where do we actually want to set/store this?
        mask = create_mask_dict(self.model.masked_latitudes,
                                self.model.masked_longitudes,
                                self.parameters['nlat'],
                                self.parameters['nlon'])
                                                            
        pruned_vtn, pruned_lambda = self.get_pruned_vtn_lambda(predict_modes)
        
        # generate KMD-based predictions
        preds = make_predictions(
            pred_horizon,
            pruned_vtn,
            pruned_lambda)
        preds[preds <= self.model.concentration_to_consider_zero] = 0
        preds[preds > 100.0] = 100.0

        np.savetxt(
            os.path.join(output_directory, 'KMD_predictions.csv'),
            preds[:self.model.number_sea_ice_observables, :],
            delimiter=",")
        n_pred = preds.shape[1]

        # reshape predictions to original longitude-latitude coordinates
        full_preds = np.zeros(
            (self.parameters['nlon'], self.parameters['nlat'], n_pred))
        for ipred in range(n_pred):
            full_preds[mask['maskedLon'], mask['maskedLat'], ipred] =\
                preds[:self.model.number_sea_ice_observables, ipred]

        dates = date_array_for_range(datetime.datetime.strptime(str(start_date),'%Y%m%d'),
                                     datetime.datetime.strptime(str(end_date), '%Y%m%d'))

        prediction_data = ClimateData(full_preds,
                                      dates,
                                      self.model.latitudes,
                                      self.model.longitudes,
                                      description=self.model.observable_name)
        #only needed if using CESM1 grid (leave off for polar stereographic)
        #prediction_data.fill_arctic_top_hole()
        return(prediction_data)
    
    def plot_snapshots(self, predictions: ClimateData,
                       mask: dict,
                       end_date: int,
                       output_directory: str):
        
        #TODO: return predictions
        
        # visualize predictios
        prediction_output_directory = os.path.join(
            output_directory, "predictions")
        os.makedirs(prediction_output_directory, exist_ok=True)
        plotPredictions(
            self.parameters,
            mask,
            predictions.latitudes,
            predictions.longitudes,
            predictions.data,
            prediction_output_directory)
    
    def select_modes_for_prediction(self, n_max: int = -1)->List[int]:
        #not 100% sure how/where I want this to live
        if n_max is not None and n_max == -1:
            predict_modes = list(range(0, self.model.modes[0].shape[2]))
        elif n_max > self.model[0].shape[2]:
            predict_modes = list(range(0, self.model.modes[0].shape[2]))
        else:
            predict_modes = list(range(0, n_max))
        return(predict_modes)
              
    def get_pruned_vtn_lambda(self,predict_modes):        
        #TODO: this might make more sense as a member of the koopman model itself.
        Koopman = {}
        Koopman['Vtn_real'] = self.model.vtn_real
        Koopman['Vtn_imag'] = self.model.vtn_imaginary
        Koopman['relambda'] = self.model.relambda
        Koopman['imlambda'] = self.model.imlambda
        Koopman['mode_imp'] = self.model.mode_imp

        # trim Koopman modes and eigenvalues based on predictModes
        Koopman['Vtn_real'] =\
            Koopman['Vtn_real'][:, predict_modes]
        Koopman['Vtn_imag'] =\
            Koopman['Vtn_imag'][:, predict_modes]
        Koopman['relambda'] =\
            Koopman['relambda'][predict_modes]
        Koopman['imlambda'] =\
            Koopman['imlambda'][predict_modes]
        Koopman['mode_imp'] =\
            Koopman['mode_imp'][predict_modes]

        # complex valued Koopman modes
        Koopman['Vtn'] = Koopman['Vtn_real'] + Koopman['Vtn_imag']*1j
        # complex valued Koopman eigenvalues
        Koopman['Lambda'] = Koopman['relambda'] + Koopman['imlambda']*1j

        m = self.model.training_data.shape[1]

        # set weights for spectral reconstruction
        W = np.ones((m, 1))

        # save coefficients of expansion and
        # scaled Koopman modes for predictions
        Koopman['coefficient'], Koopman['Vtn'] = \
            SLS_Spectral_Reconstruct_W_NE(
                Koopman['Vtn'],
                Koopman['Lambda'],
                self.model.training_data.astype(complex),
                W)
        return(Koopman['Vtn'], Koopman['Lambda'])
