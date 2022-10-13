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

import yaml
from haiku.prediction.predictor import Predictor
from haiku.training.koopman_model import KoopmanModel

if __name__ == "__main__":
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

    # load config file into parameters dictionary
    with open(args.config, 'r') as f:
        config_data = yaml.safe_load(f)

    # deserialize koopman model
    with open(args.model, 'rb') as f:
        koopman_model: KoopmanModel = pickle.load(f)

    # predict
    predictor = Predictor(config_data)

    # determine prediction modes to use
    n = args.upper_mode_bound
    if n is not None and n == -1:
        predict_modes = list(range(0, koopman_model.modes[0].shape[2]))
    else:
        predict_modes = list(range(0, n))

    predictor.predict(
        koopman_model,
        args.end_date,
        predict_modes,
        args.output_directory)
