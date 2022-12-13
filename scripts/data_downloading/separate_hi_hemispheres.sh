#!/bin/bash

# Useful for HI files which are downloaded with separate north and south datasets
# Expects filenames like "b.e11.BRCP85C5CNBDRD.f09_g16.001.cice.h.hi_nh.200601-208012.nc"
# For this example, the subdirectory will be "hi_nh"
PATTERN="\.hi_(\w{2})\."

DATA_DIRECTORY=$1

python3 split_on_pattern.py $DATA_DIRECTORY -r --pattern $PATTERN
