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
import gc

import yaml
import numpy as np
from haiku.training.trainer import Trainer
from haiku.climate_data.smear import Smearer
from haiku.climate_data.climate import load_data
from haiku.climate_data.climate_data import ClimateData
import matplotlib.pyplot as plt


def configure_logging(logging_level: Union[str, int], filename: str):
    """Create and format the global logging logger with console/file logging.

    logging_level can be: CRITICAL, ERROR, WARNING, INFO, or DEBUG
    Note: Can also use logging enums for this value (e.g. logging.INFO)
    """
    # grab the global logging logger and set the logging level
    logger = logging.getLogger()

    # ensure any other handlers are removed
    for handler in logger.handlers:
        logger.removeHandler(handler)

    logger.setLevel(logging_level)

    # create formatter to be used by the global logger
    formatter = logging.Formatter(
        '[%(asctime)s.%(msecs)03d] - %(filename)s:'
        '%(funcName)s:%(levelname)s - %(message)s',
        "%Y-%m-%d %H:%M:%S")

    # we create a stream handler for console logging
    # and a file handler for file logging
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)

    fh = logging.FileHandler(filename)
    fh.setFormatter(formatter)

    # finally, we add the console logging handler to the logger
    logger.addHandler(ch)
    logger.addHandler(fh)

if __name__ == "__main__":
    # cli args
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Path to configuration file")
    args = parser.parse_args()

    # verify
    assert os.path.exists(args.config), \
        f"{args.config} not found!"

    tracemalloc.start()
    # load config file into parameters dictionary
    with open(args.config, 'r') as f:
        config_data = yaml.safe_load(f)

    initial_config_data = copy.deepcopy(config_data)

    # TODO: Get these from configuration file
    configure_logging(logging.INFO, "haiku-training.log")

    # train koopman model given unaltered training data
    trainer = Trainer(config_data)
    trainer.train()

    #load copy of unaltered training data:
    original_climate_data = []
    print("a",trainer.parameters["serialized_data_map"])
    for serialized_data_filepath in trainer.parameters["serialized_data_map"].values():
        try:
            original_climate_data.append(ClimateData.deserialize(serialized_data_filepath))
        except FileNotFoundError as err:
            print(err)
            print("skipping this dataset")
    print("b",trainer.parameters["serialized_data_map"])
    smearer = Smearer()

    del trainer

    for i in range(config_data["n_robustness_samples"]):
        parameters = copy.deepcopy(initial_config_data)
        for key, original_data in zip(initial_config_data["serialized_data_map"],original_climate_data):
            #original_data is untouched and a smeared deep copy is returned.
            smeared_climate_data = Smearer.smear_climate_data_gaussian(original_data, parameters["smearing_factor"])
            if False:
                plt.plot(np.average(smeared_climate_data.data,(1,2)),color='k',label="smeared "+key)
                plt.plot(np.average(original_data.data,axis=(1,2)),color='green',label="measured "+key)
                plt.title(key)
                plt.legend()
                plt.show()
            n_serialized_data_path = "/local/planer_results_gallery/temp_robustness_data/"+str(i)+"_"+key+".pkl"
            smeared_climate_data.serialize(filepath=n_serialized_data_path)
            del smeared_climate_data
            #FIXME: find a better way to set this
            parameters["serialized_data_map"][key] = n_serialized_data_path


        #FIXME: Still some kind of memory leak where not everything set in trainer.train() is deallocated
        #adds colorbar label, empty clims array, nlat, nlon, nummaskobsv, numvars, NSIDCvariableName, numSeaIceObsv
        #TODO: we should save out the exact config values used during the model training somehow (make a config file that is saved alongside the koopman_model.pkl file.  And the training should run successfully start to finish if all those values are defined (instead of here where it fails due to serialized_data_map being overwritten.
        trainer = Trainer(parameters)
        print("c",trainer.parameters["serialized_data_map"])
        #FIXME: find a better way to set this
        trainer.output_directory = "/local/planer_results_gallery/temp_robustness_models/test_NSIDC_"+str(i)
        trainer.train()
        del trainer, parameters
        gc.collect()

