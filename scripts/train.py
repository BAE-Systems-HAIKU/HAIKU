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

import yaml
from haiku.training.trainer import Trainer

if __name__ == "__main__":
    # cli args
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Path to configuration file")
    args = parser.parse_args()

    # verify
    assert os.path.exists(args.config), \
        f"{args.config} not found!"

    # load config file into parameters dictionary
    with open(args.config, 'r') as f:
        config_data = yaml.safe_load(f)

    # train
    trainer = Trainer(config_data)
    trainer.train()
