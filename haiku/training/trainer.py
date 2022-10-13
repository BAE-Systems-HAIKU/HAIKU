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
import pickle
from datetime import datetime
from typing import Any, Dict

import numpy as np
import numpy.typing as npt
from haiku.climate_data.climate import (calcMask, loadCESM1Data, loadErrData,
                                        loadNSIDCData)
from haiku.climate_data.interpolation import Interpolator
from haiku.plotting.plotFunctions import plotEigenmodes, plotEigenvalues
from haiku.training.dmdR4Haiku import SLS_Spectral_Reconstruct_W_NE
from haiku.training.koopman_model import KoopmanModel
from haiku.training.training_utilities import (compute_KMD,
                                               identify_annual_variation_modes,
                                               identify_mean_modes)


class Trainer:

    def __init__(self, parameters: Dict[str, Any]):
        """Construct a Trainer object."""
        # TODO: Extract parameters into a configuration object
        self.parameters = parameters

    def train(self):
        """Initiate training processes with loaded config."""
        print(f"Loaded self.parameters: {self.parameters}")

        print('\nKoopman Mode Decomposition and Analysis on Sea Ice Data')
        print('Training on {} data'.format(self.parameters['dataFreq']))
        print('Time Interval: {} to {}'.format(
            self.parameters['start_dateInt'], self.parameters['stop_dateInt']))

        CESM1variableNameVals = self.parameters["CESM1variableNameVals"]
        NSIDCvariableNameVals = self.parameters["NSIDCvariableNameVals"]

        if np.isscalar(CESM1variableNameVals[0]):
            CESM1variableNameVals = [CESM1variableNameVals]
        if np.isscalar(NSIDCvariableNameVals[0]):
            NSIDCvariableNameVals = [NSIDCvariableNameVals]

        # loop through different experiments/datasets/files and perform KMD
        for ihemi in range(len(self.parameters["hemiVals"])):
            self.parameters['hemisphere'] = self.parameters["hemiVals"][ihemi]

            for idataSource in range(len(self.parameters["dataSourceVals"])):
                self.parameters['dataSource'] =\
                    self.parameters["dataSourceVals"][idataSource]
                if self.parameters['dataSource'] == 'CESM1':
                    # set CESM1 variable name list
                    varStrVals = CESM1variableNameVals
                    curRunTypeVals = self.parameters["runTypeVals"]
                elif self.parameters['dataSource'] == 'NSIDC':
                    # set NSIDC variable name list
                    varStrVals = NSIDCvariableNameVals
                    # NSIDC data does not have a run type
                    curRunTypeVals = ['']
                elif self.parameters['dataSource'] == 'err':
                    # This is a bad way to do it, but for 'err' runs
                    # CESM1variableNameVals and NSIDCvariableNameVals
                    # should each be length one
                    varStrVals = [1]
                    # relevant for the CESM1 data
                    curRunTypeVals = self.parameters["runTypeVals"]

                for ivarName in range(len(varStrVals)):
                    for irun in range(len(curRunTypeVals)):
                        if self.parameters['dataSource'] == 'CESM1':
                            self.parameters['CESM1variableNameVals'] =\
                                CESM1variableNameVals[ivarName]
                            numCESM1vals = len(
                                self.parameters['CESM1variableNameVals'])
                            numNSIDCvals = 0
                            self.parameters['runType'] = curRunTypeVals[irun]
                            # string with comma separated CESM1 variables names
                            varNamesStr = ','.join(
                                self.parameters['CESM1variableNameVals'])
                            print('\nRunning {} {}: {} ({})'.format(
                                self.parameters['dataSource'],
                                self.parameters['hemisphere'],
                                varStrVals[ivarName],
                                self.parameters['runType']))

                        elif self.parameters['dataSource'] == 'NSIDC':
                            self.parameters['NSIDCvariableNameVals'] =\
                                NSIDCvariableNameVals[ivarName]
                            numNSIDCvals = len(
                                self.parameters['NSIDCvariableNameVals'])
                            numCESM1vals = 0
                            # This will be '' since NSIDC
                            self.parameters['runType'] = curRunTypeVals
                            # string with comma separated NSIDC variables names
                            varNamesStr = ','.join(
                                self.parameters['NSIDCvariableNameVals'])
                            print('\nRunning {} {}: {}'.format(
                                self.parameters['dataSource'],
                                self.parameters['hemisphere'],
                                varStrVals[ivarName]))

                        elif self.parameters['dataSource'] == 'err':
                            self.parameters['CESM1variableNameVals'] =\
                                CESM1variableNameVals[ivarName]
                            self.parameters['NSIDCvariableNameVals'] =\
                                NSIDCvariableNameVals[ivarName]
                            numCESM1vals = len(
                                self.parameters['CESM1variableNameVals'])
                            numNSIDCvals = len(
                                self.parameters['NSIDCvariableNameVals'])
                            self.parameters['runType'] = curRunTypeVals[irun]
                            # string with comma separated CESM1 variables names
                            CESM1varNamesStr = ','.join(
                                self.parameters['CESM1variableNameVals'])
                            # string with comma separated NSIDC variables names
                            NSIDCvarNamesStr = ','.join(
                                self.parameters['NSIDCvariableNameVals'])
                            varNamesStr =\
                                CESM1varNamesStr + ' and ' + NSIDCvarNamesStr
                            print('\nRunning {} {}: {} ({}) and {}'.format(
                                self.parameters['dataSource'],
                                self.parameters['hemisphere'],
                                CESM1variableNameVals[ivarName],
                                self.parameters['runType'],
                                NSIDCvariableNameVals[ivarName]))

                        # set up file names if '40memberLargeEnsemble'
                        # is the chosen run type
                        if (self.parameters['runType'] ==
                                '40memberLargeEnsemble'):
                            fileNums = self.parameters['ensembleNums']
                            if np.isscalar(fileNums):
                                fileNums = [fileNums]
                        else:
                            fileNums = [0]

                        for ii in range(len(fileNums)):
                            print(
                                'Running File {} / {}'.format(
                                    ii + 1, len(fileNums)))
                            fileNum = fileNums[ii]

                            # load data
                            self.parameters['colorbar_label'] = []
                            self.parameters['clims'] = []
                            self.parameters['nlat'] = []
                            self.parameters['nlon'] = []
                            self.parameters['numMaskObsv'] = []
                            self.parameters['numVars'] = max(
                                numCESM1vals, numNSIDCvals)

                            # initialize mask dictionary
                            mask = {}
                            mask['maskedInds'] = []
                            mask['maskedLat'] = []
                            mask['maskedLon'] = []
                            mask['alphaMask'] = []

                            # loop through subvariables under ivarName
                            for isubvar in range(self.parameters['numVars']):
                                if self.parameters['dataSource'] == 'CESM1':
                                    self.parameters['CESM1variableName'] =\
                                        self.parameters[
                                            'CESM1variableNameVals'][isubvar]
                                    print('Running {} : {}'.format(
                                        self.parameters['dataSource'],
                                        self.parameters['CESM1variableName']))

                                    # Check whether the bulk version of the
                                    # current variable was specified
                                    self.parameters['isBulk'] = False
                                    if (self.parameters
                                            ['CESM1variableName'][:4]
                                            == 'bulk'):
                                        self.parameters['isBulk'] = True

                                    # data has dimensions nTime-by-nLat-by-nLon
                                    if (self.parameters
                                            ['CESM1variableName']
                                            [-len('TEMP'):]
                                            != 'TEMP'
                                            and
                                            self.parameters
                                            ['CESM1variableName'][-len('SST'):]
                                            != 'SST'
                                            and
                                            self.parameters
                                            ['CESM1variableName'][-len('HI'):]
                                            != 'HI'):
                                        # TODO: Just use climate data object
                                        # rather than extracting fields
                                        climate_data = loadCESM1Data(
                                                fileNum, self.parameters)
                                        data = climate_data.data
                                        dateInt = climate_data.date_int
                                        lat = climate_data.latitudes
                                        lon = climate_data.longitudes
                                        timeDays = climate_data.time_days
                                        descStr = climate_data.description
                                    else:
                                        # TODO: Just use climate data object
                                        # rather than extracting fields
                                        climate_data = loadCESM1Data(
                                                fileNum,
                                                self.parameters,
                                                self.parameters['seaIceLat'],
                                                self.parameters['seaIceLon'])
                                        data = climate_data.data
                                        dateInt = climate_data.date_int
                                        lat = climate_data.latitudes
                                        lon = climate_data.longitudes
                                        timeDays = climate_data.time_days
                                        descStr = climate_data.description

                                elif self.parameters['dataSource'] == 'NSIDC':
                                    self.parameters['NSIDCvariableName'] =\
                                        self.parameters[
                                            'NSIDCvariableNameVals'][isubvar]
                                    print('Running {} : {}'.format(
                                        self.parameters['dataSource'],
                                        self.parameters['NSIDCvariableName']))

                                    # Check whether the bulk version of the
                                    # current variable was specified
                                    self.parameters['isBulk'] = False
                                    if self.parameters[
                                            'NSIDCvariableName'][:4] == 'bulk':
                                        self.parameters['isBulk'] = True

                                    # data has dimensions nTime-by-nLat-by-nLon
                                    if (self.parameters['NSIDCvariableName']
                                            [-len('merged'):] == 'merged'):
                                        # TODO: Just use climate data object
                                        # rather than extracting fields
                                        climate_data = loadNSIDCData(
                                                self.parameters)
                                        data = climate_data.data
                                        dateInt = climate_data.date_int
                                        lat = climate_data.latitudes
                                        lon = climate_data.longitudes
                                        timeDays = climate_data.time_days
                                        descStr = climate_data.description
                                    else:
                                        # for sst and t2m data
                                        # (from ORAS5 and ERA5),
                                        # provide the NSIDC lat and lon arrays
                                        # to interpolate the sst or t2m data
                                        # TODO: Just use climate data object
                                        # rather than extracting fields
                                        climate_data = loadNSIDCData(
                                                self.parameters,
                                                self.parameters['seaIceLat'],
                                                self.parameters['seaIceLon'])
                                        data = climate_data.data
                                        dateInt = climate_data.date_int
                                        lat = climate_data.latitudes
                                        lon = climate_data.longitudes
                                        timeDays = climate_data.time_days
                                        descStr = climate_data.description

                                elif self.parameters['dataSource'] == 'err':
                                    # Check whether the bulk version
                                    # of the current variable was specified
                                    self.parameters['CESM1variableName'] =\
                                        self.parameters[
                                            'CESM1variableNameVals'][isubvar]
                                    self.parameters['NSIDCvariableName'] =\
                                        self.parameters[
                                            'NSIDCvariableNameVals'][isubvar]
                                    print('Running {}'.format(
                                        self.parameters['dataSource']))

                                    self.parameters['isBulk'] = False
                                    if self.parameters[
                                            'NSIDCvariableName'][:4] == 'bulk':
                                        self.parameters['isBulk'] = True
                                    # TODO: Just use climate data object
                                    # rather than extracting fields
                                    climate_data = loadErrData(
                                            fileNum, self.parameters)
                                    data = climate_data.data
                                    dateInt = climate_data.date_int
                                    lat = climate_data.latitudes
                                    lon = climate_data.longitudes
                                    timeDays = climate_data.time_days
                                    descStr = climate_data.description


                                # calculate mask indices -
                                # returns 2D coordinates
                                # [latitude_list, longitude_list]
                                # of DESIRED points
                                if (self.parameters['dataSource'] == 'CESM1' and
                                        (self.parameters['CESM1variableName'] == 'T'
                                        or
                                        self.parameters['CESM1variableName'] == 'TS'
                                        or
                                        self.parameters['CESM1variableName'] == 'flux'
                                        or
                                        self.parameters['CESM1variableName'] ==
                                            'forcing')):
                                    # for the atmosphere return all indices
                                    maskedInds = calcMask(self.parameters, -1e25, 1e25)

                                elif (self.parameters['dataSource'] == 'NSIDC' and
                                        (self.parameters['NSIDCvariableName'] == 't2m'
                                        or
                                        self.parameters['NSIDCvariableName']
                                            == 'flux')):
                                    # for the atmosphere return all indices
                                    maskedInds = calcMask(self.parameters, -1e25, 1e25)

                                else:
                                    # for all the rest,
                                    # use ice fraction coordinates
                                    maskedInds = calcMask(self.parameters, 0.0, 100.0)
                                mask['maskedInds'].append(maskedInds)

                                # these variables have already been interpolated
                                if (self.parameters['dataSource'] == 'CESM1' and
                                        (self.parameters['CESM1variableName'] == 'SST'
                                        or
                                        self.parameters['CESM1variableName'] == 'TEMP'
                                        or
                                        self.parameters['CESM1variableName'] == 'HI')):
                                    maskedLat = mask['maskedLat'][0]
                                    maskedLon = mask['maskedLon'][0]
                                else:
                                    if not self.parameters['isBulk']:
                                        if ((self.parameters['dataSource'] == 'CESM1'
                                            and
                                            (self.parameters['CESM1variableName']
                                                == 'ICEFRAC')) or
                                            (self.parameters['dataSource'] == 'NSIDC'
                                            and
                                            (self.parameters['NSIDCvariableName']
                                                == 'merged')) or
                                            (self.parameters['dataSource'] == 'err'
                                            and
                                            (self.parameters['CESM1variableName']
                                                == 'ICEFRAC')) or
                                            (self.parameters['dataSource'] == 'err'
                                            and
                                            (self.parameters['NSIDCvariableName']
                                                == 'merged'))):
                                            load_lat_num = data.shape[1]

                                    if (self.parameters['dataSource'] == 'NSIDC' or
                                            (self.parameters['dataSource'] == 'CESM1'
                                            and
                                            self.parameters['CESM1variableName']
                                                != 'forcing')):
                                        # select spatial region to use
                                        if self.parameters['hemisphere'] == 'N':
                                            inds = np.argwhere(
                                                lat >= self.parameters[
                                                    'north_lat_bound'])

                                            # flatten list of lists
                                            inds_list = [ind[0] for ind in inds]
                                            lat = lat[inds_list]
                                            data = data[:, inds_list, :]

                                            # convert latitude resolution
                                            # to degrees (-90 to 90)
                                            # latitude
                                            maskedLat = mask['maskedInds'][
                                                isubvar][0][(
                                                    -90 + 180/load_lat_num *
                                                    mask['maskedInds'][isubvar][0])
                                                            >= self.parameters[
                                                                'north_lat_bound']]
                                            # shift for plotting
                                            maskedLat = \
                                                maskedLat - load_lat_num + \
                                                inds.shape[0]

                                            # longitude
                                            maskedLon = mask['maskedInds'][
                                                isubvar][1][(
                                                    -90 + 180/load_lat_num *
                                                    mask['maskedInds'][isubvar][0])
                                                            >= self.parameters[
                                                                'north_lat_bound']]

                                        elif self.parameters['hemisphere'] == 'S':
                                            inds = np.argwhere(
                                                lat <= self.parameters[
                                                    'south_lat_bound'])

                                            # flatten list of lists
                                            inds_list = [ind[0] for ind in inds]
                                            lat = lat[inds_list]
                                            data = data[:, inds_list, :]

                                            # convert latitude resolution
                                            # to degrees (-90 to 90)
                                            # latitude
                                            maskedLat = mask['maskedInds'][
                                                isubvar][0][(
                                                    -90 + 180/load_lat_num *
                                                    mask['maskedInds'][isubvar][0])
                                                            <= self.parameters[
                                                                'south_lat_bound']]
                                            # longitude
                                            maskedLon = mask['maskedInds'][
                                                isubvar][1][(
                                                    -90 + 180/load_lat_num *
                                                    mask['maskedInds'][isubvar][0])
                                                            <= self.parameters[
                                                                'south_lat_bound']]

                                if ((not self.parameters['isBulk']) and
                                        (self.parameters['dataSource'] == 'NSIDC' or
                                        (self.parameters['dataSource'] == 'CESM1' and
                                        self.parameters['CESM1variableName']
                                            != 'forcing'))):
                                    self.parameters['nlat'].append(data.shape[1])
                                    self.parameters['nlon'].append(data.shape[2])

                                    # apply mask for plotting
                                    alphaMask = np.zeros(
                                        (data.shape[2], data.shape[1]))
                                    alphaMask[maskedLon, maskedLat] = 1.0
                                    mask['alphaMask'].append(alphaMask)

                                    # only use DESIRED points for performing KMD
                                    # data is vectorized, columns are different
                                    # snapshots and rows are different observables
                                    data = np.transpose(
                                        data[:, maskedLat, maskedLon])

                                if (self.parameters['dataSource'] == 'CESM1' and
                                        self.parameters['CESM1variableName']
                                        == 'forcing'):
                                    # only use DESIRED points for performing KMD
                                    # data is vectorized, columns are different
                                    # snapshots and rows are different observables
                                    data = np.transpose(data[:, 0:1, 0])

                                if self.parameters['isBulk']:
                                    data = np.mean(data, axis=0)
                                    data = np.expand_dims(data, axis=0)

                                if (self.parameters['isBulk'] or
                                        (self.parameters['dataSource'] == 'CESM1'
                                        and self.parameters['CESM1variableName']
                                        == 'forcing')):
                                    self.parameters['nlat'].append(1)
                                    self.parameters['nlon'].append(1)
                                    mask['maskedLat'].append(0)
                                    mask['maskedLon'].append(0)
                                else:
                                    mask['maskedLat'].append(maskedLat)
                                    mask['maskedLon'].append(maskedLon)

                                self.parameters['numMaskObsv'].append(data.shape[0])

                                # save self.parameters for sea ice
                                if not self.parameters['isBulk']:
                                    if ((self.parameters['dataSource'] == 'CESM1' and
                                        (self.parameters['CESM1variableName']
                                            == 'ICEFRAC')) or
                                        (self.parameters['dataSource'] == 'NSIDC' and
                                        (self.parameters['NSIDCvariableName']
                                            == 'merged')) or
                                        (self.parameters['dataSource'] == 'err' and
                                        (self.parameters['CESM1variableName']
                                            == 'ICEFRAC')) or
                                        (self.parameters['dataSource'] == 'err' and
                                        (self.parameters['NSIDCvariableName']
                                            == 'merged'))):
                                        self.parameters['numSeaIceObsv'] = data.shape[0]
                                        self.parameters['seaIceLat'] = lat
                                        self.parameters['seaIceLon'] = lon

                                # concatenate data matrices and create array of
                                # start and stop indices
                                # of each difference observable
                                if isubvar == 0:
                                    dataCat = data
                                    dateIntSet = dateInt
                                    indsStartDataSet = np.array([0])
                                    indsStopDataSet = np.array([np.shape(data)[0]])
                                    descStrs = [descStr]
                                    lats = [lat]
                                    lons = [lon]
                                else:
                                    dataCat = np.concatenate(
                                        (dataCat, data), axis=0)
                                    dateIntSet = np.append(
                                        dateIntSet, dateInt)
                                    indsStartDataSet = np.append(
                                        indsStartDataSet, indsStopDataSet[-1])
                                    indsStopDataSet = np.append(
                                        indsStopDataSet, indsStopDataSet[-1] +
                                        np.shape(data)[0])
                                    descStrs.append(descStr)
                                    lats.append(lat)
                                    lons.append(lon)

                                # Add colorbar label and clims for each variable
                                if (self.parameters['dataSource'] == 'CESM1' or
                                        self.parameters['dataSource'] == 'err'):
                                    if (self.parameters['CESM1variableName']
                                            [-len('ICEFRAC'):] == 'ICEFRAC'):
                                        self.parameters['colorbar_label'].append(
                                            'concentration')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['CESM1variableName']
                                            [-len('TS'):] == 'TS'):
                                        self.parameters['colorbar_label'].append(
                                            'temperature')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['CESM1variableName']
                                            [-len('T'):] == 'T'):
                                        self.parameters['colorbar_label'].append(
                                            'temperature')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['CESM1variableName']
                                            [-len('SST'):] == 'SST'):
                                        self.parameters['colorbar_label'].append(
                                            'temperature')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['CESM1variableName']
                                            [-len('TEMP'):] == 'TEMP'):
                                        self.parameters['colorbar_label'].append(
                                            'temperature')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['CESM1variableName']
                                            [-len('TEMP'):] == 'HI'):
                                        self.parameters['colorbar_label'].append(
                                            'thickness')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['CESM1variableName']
                                            [-len('flux'):] == 'flux'):
                                        self.parameters['colorbar_label'].append('flux')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['CESM1variableName']
                                            [-len('forcing'):] == 'forcing'):
                                        self.parameters['colorbar_label'].append(
                                            'forcing')
                                        self.parameters['clims'].append([])

                                elif self.parameters['dataSource'] == 'NSIDC':
                                    if (self.parameters['NSIDCvariableName']
                                            [-len('merged'):] == 'merged'):
                                        self.parameters['colorbar_label'].append(
                                            'concentration (percent)')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['NSIDCvariableName']
                                            [-len('sst'):] == 'sst'):
                                        self.parameters['colorbar_label'].append(
                                            'temperature')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['NSIDCvariableName']
                                            [-len('t2m'):] == 't2m'):
                                        self.parameters['colorbar_label'].append(
                                            'temperature')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['NSIDCvariableName']
                                            [-len('flux'):] == 'flux'):
                                        self.parameters['colorbar_label'].append('flux')
                                        self.parameters['clims'].append([])

                                    elif (self.parameters['NSIDCvariableName']
                                            [-len('thickness'):] == 'thickness'):
                                        self.parameters['colorbar_label'].append(
                                            'thickness')
                                        self.parameters['clims'].append([])

                            data = dataCat

                            print("Concatenated data shape: {}".format(data.shape))
                            mask['indsStartDataSet'] = indsStartDataSet
                            mask['indsStopDataSet'] = indsStopDataSet
                            if fileNum == 0:
                                descStrAll = varNamesStr
                            else:
                                descStrAll = 'Ensemble {}: {}'.format(
                                    fileNum, varNamesStr)

                            # create directory to save results
                            if self.parameters['runType'] == '40memberLargeEnsemble':
                                ensembleStr = '{:03d}'.format(fileNum)
                            else:
                                ensembleStr = ''

                            # setup output directory name
                            if self.parameters['dataSource'] == 'CESM1':
                                all_var_names = '+'.join(self.parameters['CESM1variableNameVals'])
                                outputStr = '{}-{}-{}_{}_{}_{}_{}to{}'.format(
                                    self.parameters['dataSource'],
                                    self.parameters['runType'],
                                    ensembleStr,
                                    all_var_names,
                                    self.parameters['hemisphere'],
                                    self.parameters['dataFreq'],
                                    self.parameters['start_dateInt'],
                                    self.parameters['stop_dateInt'])

                            elif self.parameters['dataSource'] == 'NSIDC':
                                all_var_names = '+'.join(self.parameters['NSIDCvariableNameVals'])
                                outputStr = '{}_{}_{}_{}_{}to{}'.format(
                                    self.parameters['dataSource'],
                                    all_var_names,
                                    self.parameters['hemisphere'],
                                    self.parameters['dataFreq'],
                                    self.parameters['start_dateInt'],
                                    self.parameters['stop_dateInt'])

                            elif self.parameters['dataSource'] == 'err':
                                outputStr = '{}_{}{}_and_{}_{}_{}_{}to{}'.format(
                                    self.parameters['dataSource'],
                                    '+'.join(self.parameters['CESM1variableNameVals']),
                                    ensembleStr,
                                    '+'.join(self.parameters['NSIDCvariableNameVals']),
                                    self.parameters['hemisphere'],
                                    self.parameters['dataFreq'],
                                    self.parameters['start_dateInt'],
                                    self.parameters['stop_dateInt'])
                            outputDir = os.path.join(
                                self.parameters['baseDir'],
                                self.parameters['resultsDirName'],
                                outputStr)
                            os.makedirs(outputDir, exist_ok=True)

                            # do KMD on selected data
                            Lambda, Vtn, Kefun, mode_imp = compute_KMD(
                                self.parameters["sorting"],
                                data)

                            # scaled Koopman modes
                            projCoef, Vtn = SLS_Spectral_Reconstruct_W_NE(
                                Vtn,
                                Lambda,
                                data.astype(complex),
                                np.ones((data.shape[1], 1)))

                            # convert to real and imaginary parts
                            dt = 1.0
                            relambda = Lambda.real
                            imlambda = Lambda.imag
                            omega = np.log(Lambda)/dt
                            reomega = omega.real*dt
                            imomega = (omega.imag)/(2*np.pi)

                            # create variable decoupled list of modes
                            # (last dimension is number of different variables)
                            modes = []
                            obsv_count = 0
                            for ivar in range(self.parameters['numVars']):
                                varmodes = np.zeros(
                                    (self.parameters['nlat'][ivar],
                                    self.parameters['nlon'][ivar],
                                    data.shape[1]-1))

                                for imode in range(Vtn.shape[1]):
                                    varmodes[
                                        mask['maskedLat'][ivar],
                                        mask['maskedLon'][ivar],
                                        imode] = \
                                            np.absolute(Vtn[
                                                obsv_count:obsv_count + self.parameters[
                                                    'numMaskObsv'][ivar],
                                                imode].flatten(order='C'))

                                obsv_count = \
                                    obsv_count + self.parameters['numMaskObsv'][ivar]
                                modes.append(varmodes)

                            # initialize Koopman model dictionary
                            Koopman = {
                                'relambda': relambda,
                                'imlambda': imlambda,
                                'reomega': reomega,
                                'imomega': imomega,
                                'modes': modes,
                                'mode_imp': mode_imp,
                                'Kefun': Kefun
                            }

                            koopman_model = KoopmanModel(
                                data,
                                Vtn.real,
                                Vtn.imag,
                                relambda,
                                imlambda,
                                reomega,
                                imomega,
                                modes,
                                mode_imp,
                                Kefun,
                                lat,
                                lon,
                                mask['maskedLat'][0],
                                mask['maskedLon'][0],
                                self.parameters["dataSource"],
                                all_var_names,
                                self.parameters["hemisphere"],
                                self.parameters["dataFreq"],
                                fileNum,
                                self.parameters['concentrationToConsiderZero'],
                                self.parameters['roundConc'],
                                self.parameters['north_lat_bound'],
                                self.parameters['south_lat_bound'],
                                self.parameters['numSeaIceObsv'],
                                self.parameters["start_dateInt"],
                                self.parameters["stop_dateInt"]
                            )

                            # serialize the koopman model to disk
                            koopman_model_filepath = os.path.join(
                                outputDir, "koopman_model.pkl")
                            with open(koopman_model_filepath, 'wb') as f:
                                pickle.dump(koopman_model, f)

                            # identify mean and annual variation modes
                            meanModeNums = identify_mean_modes(
                                koopman_model.reomega,
                                koopman_model.imomega,
                                koopman_model.mode_imp)
                            annualModeNums = identify_annual_variation_modes(
                                koopman_model.reomega,
                                koopman_model.imomega,
                                koopman_model.mode_imp)

                            # plot eigenfunctions and importance metric
                            # plotEigenfunctions(self.parameters, Koopman, outputDir)
                            # plotModeImportance(self.parameters, Koopman, outputDir)

                            # plot mean and annual variation modes
                            # summedMeanMode, summedAnnualMode = plotMeanAndAnnualModes(self.parameters, Koopman, mask, meanModeNums, annualModeNums, lats, lons, outputDir, descStrs)

                            # plot all eigenvalues and modes
                            if self.parameters['flag_plot_all_modes_and_eigenvalues']:
                                plotting_output_directory = os.path.join(
                                    outputDir,
                                    "plots"
                                )
                                os.makedirs(plotting_output_directory, exist_ok=True)
                                plotEigenvalues(
                                    self.parameters,
                                    Koopman,
                                    plotting_output_directory,
                                    descStrAll)
                                plotEigenmodes(
                                    self.parameters,
                                    Koopman,
                                    mask,
                                    lats,
                                    lons,
                                    plotting_output_directory,
                                    descStrs)

                            print('Output directory: {}'.format(outputDir))
