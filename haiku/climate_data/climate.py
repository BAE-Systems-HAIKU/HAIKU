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
from typing import Any, Dict

import numpy as np
import numpy.typing as npt
from haiku.climate_data.climate_data import ClimateData
from haiku.climate_data.loading_strategies.nsidc import SicMergedStrategy
from haiku.climate_data.mask import Mask


def load_mask(mask_filepath: str) -> npt.NDArray:
    """Load the preprocessed mask from a file."""
    mask = Mask.deserialize(mask_filepath)
    return mask.indices


def load_data(serialized_filepath: str, start_date: int,
              stop_date: int, round_concentration: bool,
              concentration_to_consider_zero: float) -> ClimateData:
    """Load Dataset as ClimateData.

    Deserializes the preprocessed data and filters it.
    """
    #FIXME: assumes climate data has shape: [ntime, nlat, nlon]
    #FIXME: would prefer if climate data object had shape: [nlon, nlat, ntime]
    # deserialize data from disk
    assert os.path.exists(serialized_filepath), \
        f"{serialized_filepath} not found"
    climate_data = ClimateData.deserialize(serialized_filepath)

    # TODO: Refactor filters
    # filter by date

    # find the index of the closest elements from the array for start/stop
    # in the difference array
    # this is in case the index is not present in the array
    #TODO: should we return an error if the requested stop date is
    #more than one time step beyond the final time in dataset?
    lowest_index = np.absolute(climate_data.date_int - start_date).argmin()
    highest_index = np.absolute(climate_data.date_int - stop_date).argmin()

    # TODO: Why is there +1 here
    #Answer: we want to be inclusive of both start and end dates
    timeInds = range(lowest_index, highest_index + 1)
    climate_data.data = climate_data.data[timeInds, :, :]
    climate_data.date_int = climate_data.date_int[timeInds]
    climate_data.dates = climate_data.dates[timeInds]
    
    # round to nearest percent
    if round_concentration:
        climate_data.data = np.around(climate_data.data)

    # sea ice concentration percentages
    # at or below this value will be set to zero
    # also check that max values are capped at 100%
    # TODO: Ensure this only happens if the datasource is not err
    if climate_data.description == "ICEFRAC" or \
            climate_data.description == SicMergedStrategy.VARIABLE_KEY:
        climate_data.data[
            climate_data.data <= concentration_to_consider_zero] = 0.0

    return climate_data


def loadErrData(fileNum: str, parameters: Dict[str, Any]) -> ClimateData:
    """Load difference of CESM1 and NSIDC Datasets into Numpy arrays."""
    # Load CESM1 data
    parameters['dataSource'] = 'CESM1'
    dataC, dateIntC, latC, lonC, timeDaysC, descStrC =\
        loadCESM1Data(fileNum, parameters)

    # Load NSIDC data
    parameters['dataSource'] = 'NSIDC'
    dataN, dateIntN, latN, lonN, timeDaysN, descStrN =\
        loadNSIDCData(parameters)

    parameters['dataSource'] = 'err'

    # Compute error
    data = dataC - dataN
    descStr = 'Error CESM1 {} - NSIDC {}'.format(
        parameters['CESM1variableName'],
        parameters['NSIDCvariableName'])

    dateInt = dateIntC
    lat = latC
    lon = lonC
    timeDays = timeDaysC

    return ClimateData(data, dateInt, lat, lon, timeDays, descStr)


def climMeanNSIDC(parameters: Dict[str, Any]) -> npt.NDArray:
    """Compute the climatological mean using the NSIDC dataset."""
    # create temporary parameter dictionary for calling loadNSIDCData
    #TODO: we shouldn't need to set temporary parameters, we can just
    #use the relevant ones from NSIDC object and then drop them 
    temp_parameters = {}
    temp_parameters['dataSource'] = parameters['dataSource']
    temp_parameters['NSIDCvariableName'] = parameters['NSIDCvariableName']
    temp_parameters['hemisphere'] = parameters['hemisphere']
    temp_parameters['start_dateInt'] = 19781101
    temp_parameters['stop_dateInt'] = 20201201
    # do not use artificial data
    temp_parameters['flag_use_missing'] = False
    # not a bulk variable
    temp_parameters['isBulk'] = False
    # round sea ice concentration
    temp_parameters['roundConc'] = True
    temp_parameters['concentrationToConsiderZero'] = 15

    # load the complete NSIDC dataset
    serialized_data_filepath =\
        parameters["serialized_data_map"][
            temp_parameters['NSIDCvariableName']]
    climate_data = load_data(
        serialized_data_filepath,
        temp_parameters["start_dateInt"],
        temp_parameters["stop_dateInt"],
        parameters["roundConc"],
        parameters["concentrationToConsiderZero"]
    )
    comp_data = climate_data.data
    dateInt = climate_data.date_int

    comp_data[comp_data <= parameters['concentrationToConsiderZero']] = 0

    # retrieve the corresponding month for each snapshot
    monthInt = dateInt // 100 % 100

    # initialize climatological mean
    climMean = np.zeros((12, comp_data.shape[1], comp_data.shape[2]))
    for imonth in range(12):
        # take the average of one month over the entire dataset
        climMean[imonth, :, :] = np.mean(
            comp_data[monthInt == (imonth+1), :, :], axis=0)

    return climMean
