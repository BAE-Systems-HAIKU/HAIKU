#!/bin/bash

# Expects filenames like "b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h0.ICEFRAC.185001-200512.nc"
# In this example, the file will be put in a subdirectory balled "001"
PATTERN="\.\d{3}\."

DATA_DIRECTORY=$1

python3 split_on_pattern.py $DATA_DIRECTORY -r --pattern $PATTERN
