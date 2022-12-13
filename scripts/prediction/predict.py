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
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

import yaml
from haiku.prediction.predictor import Predictor
from haiku.training.koopman_model import KoopmanModel
from haiku.climate_data.climate import load_data
from haiku.climate_data.climate_data import ClimateData
from haiku.utilities.mask import create_mask_dict
from haiku.plotting.plotFunctions import (plot_timeseries, plotSnapshot, plotSnapshotsSideBySide)



def arg_parse():
    # cli args
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Path to configuration file")
    parser.add_argument("model", type=str, help="Path to serialized model")
    parser.add_argument("end_date", type=int, help="YYYYMMDD (e.g., 20101201")
    parser.add_argument("output_directory", type=str,
                        help="Where to output plots")
    parser.add_argument("--upper_mode_bound", type=int, default=-1,
                        help="""Prediction will use modes 0 to this upperbound.
                        If not specified, all modes will be used.""")
    args = parser.parse_args()

    # verify
    assert os.path.exists(args.config), \
        f"{args.config} not found!"
    assert os.path.exists(args.model), \
        f"{args.model} not found!"
    os.makedirs(args.output_directory, exist_ok=True)

    return(args)


def main(args):

    # load config file into parameters dictionary
    with open(args.config, 'r') as f:
        config_data = yaml.safe_load(f)

    # deserialize koopman model
    with open(args.model, 'rb') as f:
        koopman_model: KoopmanModel = pickle.load(f)

    # predict
    predictor = Predictor(koopman_model)

    # determine prediction modes to use
    n = args.upper_mode_bound
    if n is not None and n == -1:
        predict_modes = list(range(0, koopman_model.modes[0].shape[2]))
    else:
        predict_modes = list(range(0, n))

    start_date = predictor.model.training_start_date
    prediction = predictor.predict(
        start_date,
        args.end_date,
        predict_modes,
        args.output_directory)

    mask = create_mask_dict(koopman_model.masked_latitudes,
                            koopman_model.masked_longitudes,
                            len(koopman_model.latitudes),
                            len(koopman_model.longitudes))
    NSIDC_data = load_data(
        config_data['serialized_data_map']["merged"],
        start_date,
        args.end_date,
        config_data["roundConc"],
        config_data["concentrationToConsiderZero"])

    climate_mean_data = NSIDC_data.return_climatological_mean()

    CESM_data = load_data(
        config_data['serialized_data_map']["ICEFRAC"],
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
        reference_data[label].select_hemisphere("POLAR")
        #TODO: update climate data storage so that transpose isn't needed
        reference_data[label].data = np.transpose(reference_data[label].data)
        for iday in range(reference_data[label].data.shape[2]):
            reference_data[label].data[:,:,iday] = reference_data[label].data[:,:,iday]*mask["alphaMask"]
        #set top hole in arctic to 100%
        #only needed if using CESM1 grid (leave off for polar stereographic)
        #reference_data[label].fill_arctic_top_hole()
        #mask already applied to prediction (in predict function) and reference data
    plot_timeseries(reference_data,
                    prediction,
                    args.output_directory,
                    "coverage")
    plot_timeseries(reference_data,
                    prediction,
                    args.output_directory,
                    "spatial_correlation")
    plot_timeseries(reference_data,
                    prediction,
                    args.output_directory,
                    "rmse")

    os.makedirs(args.output_directory+"predictions/", exist_ok=True)
    os.makedirs(args.output_directory+"comparisons/", exist_ok=True)
    for i,date in enumerate(prediction.date_int):
        #generate 2-d spatial data for predictions
        image_title = "KMD_prediction_"+str(date)
        fig = plt.figure()
        fig = plotSnapshot(fig,
                           prediction.data[:,:,i],
                           image_title,
                           "Sea Ice Percentage",
                           [0,100],
                           prediction.latitudes,
                           prediction.longitudes,
                           mask["alphaMask"])

        plt.savefig(args.output_directory+"predictions/"+image_title+".png", bbox_inches='tight', pad_inches=.1)
        plt.close("all")

        #generate comparisons of each timestep 2-d spatial data
        image_title = "comparison_KMD_prediction_"+str(date)
        fig = plt.figure()
        fig = plotSnapshotsSideBySide(fig,
                                     reference_data["NSIDC"].data[:,:,i],
                                     prediction.data[:,:,i],
                                     image_title,
                                     "Sea Ice Percentage",
                                     [0,100],
                                     prediction.latitudes,
                                     prediction.longitudes,
                                     mask["alphaMask"])
        plt.savefig(args.output_directory+"comparisons/"+image_title+".png", bbox_inches='tight', pad_inches=.1)
        plt.close("all")


if __name__ == "__main__":
    args = arg_parse()
    main(args)
