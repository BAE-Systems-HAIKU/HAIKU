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
from typing import Union

import yaml
from haiku.training.trainer import Trainer


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


def main(args):

    # verify
    assert os.path.exists(args.config), \
        f"{args.config} not found!"

    # load config file into parameters dictionary
    with open(args.config, 'r') as f:
        config_data = yaml.safe_load(f)

    # configure logging (put log in same directory as script)
    log_filepath =\
        os.path.join(os.path.dirname(__file__), "haiku_training.log")
    configure_logging(logging.INFO, log_filepath)

    # train
    trainer = Trainer(config_data)
    trainer.train()
    
if __name__ == "__main__":
    # cli args
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Path to configuration file")
    args = parser.parse_args()
    main(args)
