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

import argparse
import logging
import os
import copy
from typing import Union
import tracemalloc
from pathlib import Path
import numpy as np
import gc

import yaml
import pickle
from haiku.training.trainer import Trainer
from haiku.climate_data.smear import Smearer
from haiku.climate_data.climate import load_data
from haiku.climate_data.climate_data import ClimateData
from haiku.utilities.mask import create_mask_dict
import matplotlib.pyplot as plt

from haiku.prediction.predictor import Predictor
from haiku.training.koopman_model import KoopmanModel
from haiku.plotting.plotFunctions import (plotTimeSeries,plot_robustness_timeseries, plot_timeseries, plotSnapshotsSideBySide)


if __name__ == "__main__":
    # cli args                                                                                                        
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Path to configuration file")
    parser.add_argument("models_path", type=str, help="Path to directory holding subdirectories with koopman_model.pkl")
    parser.add_argument("original_model", type=str, help="path to directory holding baseline koopman model to compare with modified koopman models for robustness analysis")
    parser.add_argument("end_date", type=int, help="YYYYMMDD (e.g., 20101201")
    parser.add_argument("output_directory", type=str, help="Where to output plots")
    parser.add_argument("--upper_mode_bound", type=int, default=-1,
                        help="""Prediction will use modes 0 to this upperbound.
                        If not specified, all modes will be used.""")
    parser.add_argument("--target_dataset", type=str, default="NSIDC",
                        help="""specific the target dataset to compare with
                        defaults to NSIDC. Uses serialized_data_maps and 
                        XYZvariableNameVals from config file""")
    parser.add_argument("--additional_climate_data", type=list, default=["Climatalogical Mean"],
                         help="""list of any additional datasets to be loaded from
                         climate data pkl objects to be included on plots (will always plot sea ice)
                         defaults to empty array. Uses serialized_data_maps and 
                        XYZvariableNameVals from config file""")
                        
    args = parser.parse_args()
    

    # verify
    assert os.path.exists(args.config), \
        f"{args.config} not found!"
    assert os.path.exists(args.models_path), \
        f"{args.models_path} not found!"
    assert os.path.exists(args.original_model), \
        f"{args.original_model} not found!"
    os.makedirs(args.output_directory, exist_ok=True)

    # load config file into parameters dictionary
    with open(args.config, 'r') as f:
        config_data = yaml.safe_load(f)
        
    model_paths = list(Path(args.models_path).rglob("*.pkl"))
    predictions = []


    ######## Load model from original data for comparison ######
    with open(args.original_model, 'rb') as f:
        original_model: KoopmanModel = pickle.load(f)

    # predict
    predictor = Predictor(original_model)
    start_date = predictor.model.training_start_date
    predict_modes = predictor.select_modes_for_prediction(args.upper_mode_bound)
    original_prediction = predictor.predict(
        start_date,
        args.end_date,
	predict_modes,
	args.output_directory+"/original")

    ######## Create mask from koopman model parameters ###########
    #latitudes and longitudes really represent generic x, y coordinates in this case
    mask = create_mask_dict(original_model.masked_latitudes,
                            original_model.masked_longitudes,
			    len(original_model.latitudes),
                            len(original_model.longitudes))
    #Load the NSIDC data
    NSIDC_data = load_data(
        config_data['serialized_data_map'][config_data['NSIDCvariableName']],
        start_date,
        args.end_date,
        config_data["roundConc"],
        config_data["concentrationToConsiderZero"])

    #compute the climatological mean
    climate_mean_data = NSIDC_data.return_climatological_mean()

    #load the CESM1 data
    CESM_data = load_data(
        config_data['serialized_data_map'][config_data['CESM1variableName']],
        start_date,
        args.end_date,
        config_data["roundConc"],
        config_data["concentrationToConsiderZero"])

    reference_data = {}
    monthInt_nsidc = (NSIDC_data.date_int // 100 % 100) -1
    yearMonth_nsidc = (NSIDC_data.date_int // 100)

    for temp_data, label in zip([NSIDC_data, climate_mean_data, CESM_data],
                         ["NSIDC", "climatological mean", "CESM1"]):
        reference_data[label] = temp_data
        reference_data[label].select_hemisphere(config_data['hemisphere'])
        #TODO: update climate data storage so that transpose isn't needed
        reference_data[label].data = np.transpose(reference_data[label].data)
        for iday in range(reference_data[label].data.shape[2]):
            reference_data[label].data[:,:,iday] = reference_data[label].data[:,:,iday]*mask["alphaMask"]

        #set top hole in arctic to 100%
	#only needed if using CESM1 grid (leave off for polar stereographic)
        #reference_data[label].fill_arctic_top_hole()
        #mask already applied to prediction (in predict function) and reference data

    #plots of default model:
    plot_timeseries(reference_data,
                    original_prediction,
                    args.output_directory+"/original/",
                    "coverage")
    plot_timeseries(reference_data,
                    original_prediction,
                    args.output_directory+"/original/",
                    "spatial_correlation")
    plot_timeseries(reference_data,
                    original_prediction,
                    args.output_directory+"/original/",
                    "rmse")

    predictor.parameters = config_data
    predictor.parameters["dataSource"]=predictor.model.data_source
        
    #TODO: check that this step is unneeded:
    #should already have entries for each timestep in NSIDC dataset
    #reference_data["climatological_mean"] = reference_data["climatological_mean"][:,:, monthInt_nsidc]
        
    for i, model in enumerate(model_paths):
        # deserialize koopman model
        with open(model, 'rb') as f:
            koopman_model: KoopmanModel = pickle.load(f)

        # predict
        predictor = Predictor(koopman_model)

        # determine prediction modes to use
        predict_modes = predictor.select_modes_for_prediction(args.upper_mode_bound)
        os.makedirs(args.output_directory+"/"+str(i), exist_ok=True)
        predictions.append(predictor.predict(
            predictor.model.training_start_date,
            args.end_date,
	    predict_modes,
	    args.output_directory+"/"+str(i)))

        plot_timeseries(reference_data,
                        predictions[i],
                        args.output_directory+"/"+str(i),
                        "coverage")
        plot_timeseries(reference_data,
                        predictions[i],
                        args.output_directory+"/"+str(i),
                        "spatial_correlation")
        #RMSE doesn't make much sense since the perturbed datasets and models
        #should always be worse than the original model predictions vs original dataset
        plot_timeseries(reference_data,
                        predictions[i],
                        args.output_directory+"/"+str(i),
                        "rmse")
        #add plotsnapshotsidebyside() if you want to see each model's prediction at each time step
        #on the 2-d spatial grid compared to the unadultered training data.
        #generate comparisons of each timestep 2-d spatial data


        if config_data["flag_plot_comparisons"] == True:
            os.makedirs(args.output_directory+"comparisons/", exist_ok=True)
            for i_time, date in enumerate(predictions[i].date_int):
                image_title = "comparison_KMD_prediction_"+str(date)
                fig = plt.figure()
                fig = plotSnapshotsSideBySide(fig,
                                              reference_data["NSIDC"].data[:,:,i_time],
                                              predictions[i].data[:,:,i_time],
                                              image_title,
                                              "Sea Ice Percentage",
                                              [0,100],
                                              predictions[i].latitudes,
                                              predictions[i].longitudes,
                                              mask["alphaMask"])
                plt.savefig(args.output_directory+"comparisons/"+image_title+".png",
                            bbox_inches='tight', pad_inches=.1)
                plt.close("all")
                
                #copy the original Koopman model prediction for shape
                #compute the average of all the koopman model predictions
                #plot the average of all the koopman model predictions with robustness bands
    average_prediction = copy.deepcopy(original_prediction)
                
    average_prediction.average(predictions)
    average_prediction.compute_sigmas(predictions)
    
    plot_robustness_timeseries(reference_data,
                               average_prediction,
                               args.output_directory+"/average/",
                               "coverage")
    plot_robustness_timeseries(reference_data,
                               average_prediction,
                               args.output_directory+"/average/",
                               "spatial_correlation")
    #RMSE doesn't make much sense since the perturbed datasets and models
    #should always be worse than the original model predictions vs original dataset
    plot_robustness_timeseries(reference_data,
                               average_prediction,
                               args.output_directory+"/average/",
                               "rmse")
    
