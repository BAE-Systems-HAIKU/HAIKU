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

import logging
import os
import pickle
from typing import Any, Dict, List, Mapping

import numpy as np
import numpy.typing as npt
from haiku.climate_data.climate import load_data, load_mask, loadErrData
from haiku.plotting.plotFunctions import plotEigenmodes, plotEigenvalues
from haiku.training.dmdR4Haiku import SLS_Spectral_Reconstruct_W_NE
from haiku.training.koopman_model import KoopmanModel
from haiku.training.training_utilities import compute_KMD


class Trainer:

    def __init__(self, parameters: Dict[str, Any]):
        """Construct a Trainer object."""
        # TODO: Extract parameters into a configuration object
        self.parameters = parameters
        self.data_sources = self.parameters["data_sources"]
        self.dataset_paths = []
        
    def train(self):
        """Initiate training processes with loaded config."""
        logging.info("Starting training method...")
        logging.info("Loaded parameters: %s", self.parameters)
        logging.info("Time Interval: %s to %s",
                     self.parameters['start_dateInt'],
                     self.parameters['stop_dateInt'])

        # TODO: This is a bit odd, let's simplify at some point
        CESM1variableNameVals = self.parameters["CESM1variableNameVals"]
        NSIDCvariableNameVals = self.parameters["NSIDCvariableNameVals"]
        if np.isscalar(CESM1variableNameVals[0]):
            CESM1variableNameVals = [CESM1variableNameVals]
        if np.isscalar(NSIDCvariableNameVals[0]):
            NSIDCvariableNameVals = [NSIDCvariableNameVals]

        # perform KMD on configured data source
        for ds in self.data_sources:
            self.data_source = ds
            logging.info("Processing data source: %s", self.data_source)
            self._process_data_source(
                CESM1variableNameVals,
                NSIDCvariableNameVals)

    def _process_data_source(self, CESM1variableNameVals: List[str],
                             NSIDCvariableNameVals: List[str]):
        # determine variables and run type from data source
        if self.data_source == 'CESM1':
            # set CESM1 variable name list
            varStrVals = CESM1variableNameVals
        elif self.data_source == 'NSIDC':
            # set NSIDC variable name list
            varStrVals = NSIDCvariableNameVals
        elif self.data_source == 'err':
            # This is a bad way to do it, but for 'err' runs
            # CESM1variableNameVals and NSIDCvariableNameVals
            # should each be length one
            varStrVals = [1]

        # process each variable and run type for data source
        for ivarName in range(len(varStrVals)):
            self._process_variable(
                ivarName,
                CESM1variableNameVals,
                NSIDCvariableNameVals,
                varStrVals)

    def _process_variable(self,
                          ivarName: int,
                          CESM1variableNameVals: List[str],
                          NSIDCvariableNameVals: List[str],
                          varStrVals: List[str]):
        if self.data_source == 'CESM1':
            self.parameters['CESM1variableNameVals'] =\
                CESM1variableNameVals[ivarName]
            numCESM1vals = len(
                self.parameters['CESM1variableNameVals'])
            numNSIDCvals = 0
            # string with comma separated CESM1 variables names
            varNamesStr = ','.join(
                self.parameters['CESM1variableNameVals'])
            print('\nRunning {} {}: {}'.format(
                self.data_source,
                self.parameters['hemisphere'],
                varStrVals[ivarName]))

        elif self.data_source == 'NSIDC':
            self.parameters['NSIDCvariableNameVals'] =\
                NSIDCvariableNameVals[ivarName]
            numNSIDCvals = len(
                self.parameters['NSIDCvariableNameVals'])
            numCESM1vals = 0
            # string with comma separated NSIDC variables names
            varNamesStr = ','.join(
                self.parameters['NSIDCvariableNameVals'])
            print('\nRunning {} {}: {}'.format(
                self.data_source,
                self.parameters['hemisphere'],
                varStrVals[ivarName]))

        elif self.data_source == 'err':
            self.parameters['CESM1variableNameVals'] =\
                CESM1variableNameVals[ivarName]
            self.parameters['NSIDCvariableNameVals'] =\
                NSIDCvariableNameVals[ivarName]
            numCESM1vals = len(
                self.parameters['CESM1variableNameVals'])
            numNSIDCvals = len(
                self.parameters['NSIDCvariableNameVals'])
            # string with comma separated CESM1 variables names
            CESM1varNamesStr = ','.join(
                self.parameters['CESM1variableNameVals'])
            # string with comma separated NSIDC variables names
            NSIDCvarNamesStr = ','.join(
                self.parameters['NSIDCvariableNameVals'])
            varNamesStr =\
                CESM1varNamesStr + ' and ' + NSIDCvarNamesStr
            print('\nRunning {} {}: {} and {}'.format(
                self.data_source,
                self.parameters['hemisphere'],
                CESM1variableNameVals[ivarName],
                NSIDCvariableNameVals[ivarName]))

        # TODO: Remove filenums-- this doesn't really do anything anymore
        fileNums = [0]
        for ii in range(len(fileNums)):
            self._process_file_number(
                ii, fileNums, numCESM1vals, numNSIDCvals, varNamesStr)

    def _process_file_number(self, ii: int, fileNums: List[int],
                             numCESM1vals: int, numNSIDCvals: int,
                             varNamesStr: str):
        print('Running File {} / {}'.format(ii + 1, len(fileNums)))
        fileNum = fileNums[ii]

        # load data
        self.parameters['colorbar_label'] = []
        self.parameters['clims'] = []
        self.parameters['nlat'] = []
        self.parameters['nlon'] = []
        self.parameters['numMaskObsv'] = []
        self.parameters['numVars'] = max(numCESM1vals, numNSIDCvals)

        # initialize mask dictionary
        mask = {}
        mask['maskedInds'] = []
        mask['maskedLat'] = []
        mask['maskedLon'] = []
        mask['alphaMask'] = []

        # loop through subvariables under ivarName
        for isubvar in range(self.parameters['numVars']):
            data, dateInt, lat, lon, descStr =\
                self._process_subvariable(isubvar, fileNum, varNamesStr, mask)

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

            self._add_colorbar_labels_and_clims()

        data = dataCat

        print("Concatenated data shape: {}".format(data.shape))
        mask['indsStartDataSet'] = indsStartDataSet
        mask['indsStopDataSet'] = indsStopDataSet
        if fileNum == 0:
            descStrAll = varNamesStr
        else:
            descStrAll = 'Ensemble {}: {}'.format(
                fileNum, varNamesStr)

        outputDir, all_var_names =\
            self._determine_output_directory()
        os.makedirs(outputDir, exist_ok=True)

        self._train_save_and_plot(data, mask, lat, lon, all_var_names, fileNum,
                                  outputDir, descStrs, descStrAll, lats, lons)

    def _process_subvariable(self, isubvar: int,
                             fileNum: int, varNamesStr: str,
                             mask: Mapping[str, Any]):
            
        if self.data_source == 'CESM1':
            var_name_key = 'CESM1variableName'
            var_name_val_key = 'CESM1variableNameVals'
        elif self.data_source == 'NSIDC':
            var_name_key = 'NSIDCvariableName'
            var_name_val_key = 'NSIDCvariableNameVals'
            
        if self.data_source == 'NSIDC' or self.data_source == 'CESM1':
            self.parameters[var_name_key] =\
                self.parameters[
                    var_name_val_key][isubvar]
            print('Running {} : {}'.format(
                self.data_source,
                self.parameters[var_name_key]))

            maskedInds = self._determine_mask_indices()
            mask['maskedInds'].append(maskedInds)

            # Check whether the bulk version of the
            # current variable was specified
            self.parameters['isBulk'] = False
            if self.parameters[
                    var_name_key][:4] == 'bulk':
                self.parameters['isBulk'] = True

            # data has dimensions nTime-by-nLat-by-nLon
            # TODO: Dynamically determine data filepath
            serialized_data_filepath =\
                self.parameters["serialized_data_map"][
                    self.parameters[var_name_key]]
            climate_data = load_data(
                serialized_data_filepath,
                self.parameters["start_dateInt"],
                self.parameters["stop_dateInt"],
                self.parameters["roundConc"],
                self.parameters["concentrationToConsiderZero"]
            )
            data = climate_data.data
            dateInt = climate_data.date_int
            lat = climate_data.latitudes
            lon = climate_data.longitudes
            descStr = climate_data.description

        elif self.data_source == 'err':
            # Check whether the bulk version
            # of the current variable was specified
            self.parameters['CESM1variableName'] =\
                self.parameters[
                    'CESM1variableNameVals'][isubvar]
            self.parameters['NSIDCvariableName'] =\
                self.parameters[
                    'NSIDCvariableNameVals'][isubvar]
            print('Running {}'.format(
                self.data_source))

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
            descStr = climate_data.description

        # these variables have already been interpolated
        if (self.data_source == 'CESM1' and
                (self.parameters['CESM1variableName'] == 'SST'
                    or
                    self.parameters['CESM1variableName'] == 'TEMP'
                    or
                    self.parameters['CESM1variableName'] == 'HI')):
            maskedLat = mask['maskedLat'][0]
            maskedLon = mask['maskedLon'][0]
            load_lat_num = data.shape[1]
        else:
            if not self.parameters['isBulk']:
                load_lat_num = data.shape[1]

            if (self.data_source == 'NSIDC' or
                    (self.data_source == 'CESM1'
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
                elif self.parameters['hemisphere'] == 'POLAR':
                    # hacky way to get all inds
                    # TODO: There's probably an np way to get this
                    inds = np.argwhere(lat == lat)

                    # flatten list of lists
                    inds_list = [ind[0] for ind in inds]
                    lat = lat[inds_list]
                    data = data[:, inds_list, :]

                    maskedLat = mask['maskedInds'][isubvar][0]
                    maskedLon = mask['maskedInds'][isubvar][1]
                else:
                    raise ValueError("Invalid Hemisphere provided: %s"
                                     % self.parameters['hemisphere'])

        if ((not self.parameters['isBulk']) and
                (self.data_source == 'NSIDC' or
                    (self.data_source == 'CESM1' and
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

        if (self.data_source == 'CESM1' and
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
                (self.data_source == 'CESM1'
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
            if ((self.data_source == 'CESM1' and
                (self.parameters['CESM1variableName']
                    == 'ICEFRAC')) or
                (self.data_source == 'NSIDC' and
                (self.parameters['NSIDCvariableName']
                    == 'merged')) or
                (self.data_source == 'err' and
                (self.parameters['CESM1variableName']
                    == 'ICEFRAC')) or
                (self.data_source == 'err' and
                (self.parameters['NSIDCvariableName']
                    == 'merged'))):
                self.parameters['numSeaIceObsv'] = data.shape[0]
                self.parameters['seaIceLat'] = lat
                self.parameters['seaIceLon'] = lon

        print(serialized_data_filepath)
        print(data.shape)
        
        return data, dateInt, lat, lon, descStr

    def _train_save_and_plot(self, data: npt.NDArray, mask: Mapping,
                             lat: npt.NDArray, lon: npt.NDArray,
                             all_var_names: str, fileNum: int,
                             outputDir: str, descStrs: str, descStrAll: str,
                             lats, lons):
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
            varmodes = np.zeros((
                self.parameters['nlat'][ivar],
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
            self.data_source,
            all_var_names,
            self.parameters["hemisphere"],
            fileNum,
            "Monthly",
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
            print("parameters in place for plotting",str(self.parameters["numVars"]))
            print(self.parameters["colorbar_label"])
            plotEigenmodes(
                self.parameters,
                Koopman,
                mask,
                lats,
                lons,
                plotting_output_directory,
                descStrs)

        print('Output directory: {}'.format(outputDir))

    def _determine_mask_indices(self):
        return load_mask(self.parameters["mask_filepath"])

    def _add_colorbar_labels_and_clims(self):
        # Add colorbar label and clims for each variable
        if (self.data_source == 'CESM1' or
                self.data_source == 'err'):
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

            elif (self.parameters['CESM1variableName'][-len('TEMP'):] == 'HI' or
                    self.parameters['CESM1variableName'][-len('TEMP'):] == 'hi'):
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

            else:
                raise ValueError(
                    "No Colorbar label entry to variable: %s"
                    % self.parameters['CESM1variableName'])

        elif self.data_source == 'NSIDC':
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

            else:
                raise ValueError(
                    "No Colorbar label entry to variable: %s"
                    % self.parameters['NSIDCvariableName'])

    def _determine_output_directory(self) -> str:

        #allow custom naming rather than predefined structure
        try:
            type(self.output_directory)
            if self.data_source == "CESM1":
                all_var_names = '+'.join(self.parameters['CESM1variableNameVals'])
            elif self.data_source == "NSIDC":
                all_var_names = '+'.join(self.parameters['CESM1variableNameVals'])
            return(self.output_directory, all_var_names)
        #if self.output_directory was undefined at instantiation, follow standard structure from config file
        except AttributeError:
            pass
        
        if self.data_source == 'CESM1':
            all_var_names = '+'.join(self.parameters['CESM1variableNameVals'])
            outputStr = '{}-{}-{}_{}to{}'.format(
                self.data_source,
                all_var_names,
                self.parameters['hemisphere'],
                self.parameters['start_dateInt'],
                self.parameters['stop_dateInt'])

        elif self.data_source == 'NSIDC':
            all_var_names = '+'.join(self.parameters['NSIDCvariableNameVals'])
            outputStr = '{}_{}_{}_{}to{}'.format(
                self.data_source,
                all_var_names,
                self.parameters['hemisphere'],
                self.parameters['start_dateInt'],
                self.parameters['stop_dateInt'])

        elif self.data_source == 'err':
            outputStr = '{}_{}_and_{}_{}_{}to{}'.format(
                self.data_source,
                '+'.join(self.parameters['CESM1variableNameVals']),
                '+'.join(self.parameters['NSIDCvariableNameVals']),
                self.parameters['hemisphere'],
                self.parameters['start_dateInt'],
                self.parameters['stop_dateInt'])
        outputDir = os.path.join(
            self.parameters['output_directory'],
            outputStr)
        return outputDir, all_var_names
