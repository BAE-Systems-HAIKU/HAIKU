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

import glob
import os
from typing import Any, Dict

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from haiku.climate_data.climate_data import ClimateData
from haiku.climate_data.interpolation import Interpolator
from haiku.climate_data.readCESM1 import (readCESM1FLNT, readCESM1forcing,
                                          readCESM1FSNT, readCESM1HI,
                                          readCESM1ICEFRAC, readCESM1SST,
                                          readCESM1T, readCESM1TEMP,
                                          readCESM1TS)
from haiku.climate_data.readNSIDC import (read_ERA5_surf_flux, read_ERA5_t2m,
                                          read_ERA5_top_flux, read_ORAS5_sst,
                                          read_PSC_thickness, readNSIDCSIC)
from netCDF4 import Dataset
from scipy import interpolate


def calcMask(parameters: Dict[str, Any], min: float, max: float,
             visualize_masking: bool = False) -> npt.NDArray:
    """Compute masked and unmasked indices from earliest NSIDC data.

    Earliest NSIDC data is (1978-11).

    """
    # pull variables from parameter dictionary
    baseDir = parameters['baseDir']
    hemisphere = parameters['hemisphere']

    varName = 'goddard_merged_seaice_conc_monthly'

    #TODO: can we just use the saved grid files here?
    if hemisphere == 'N':
        fn = 'latlon_seaice_conc_monthly_nh_n07_197811_v03r01.nc'
        fnFull = os.path.join(
            baseDir,
            parameters['dataDirName'],
            'NSIDC/monthly/north_regrid/merged',
            fn)
    elif hemisphere == 'S':
        fn = 'latlon_seaice_conc_monthly_sh_n07_197811_v03r01.nc'
        fnFull = os.path.join(
            baseDir,
            parameters['dataDirName'],
            'NSIDC/monthly/south_regrid/merged',
            fn)

    # open NetCDF file
    ncfile = Dataset(fnFull)
    mapping = np.squeeze(
        np.array(ncfile.variables[varName]))
    ncfile.close()

    # obtain vectorized mask indices
    # contains all DESIRED points
    M = np.where((mapping >= min) & (mapping <= max))

    # cast to numpy array
    M = np.asarray(M)

    if visualize_masking:
        _ = plt.figure(1)
        plt.title("Unmasked Map")
        plt.imshow(np.flipud(mapping))

        # enforce masking pixel value
        mapping[M[0], M[1]] = -100

        _ = plt.figure(2)
        temp_plot = np.ma.masked_where((mapping == -100), mapping)
        cmap = plt.cm.get_cmap("viridis").copy()
        cmap.set_bad(color='red')
        plt.imshow(np.flipud(temp_plot), cmap=cmap)
        plt.title("Masked Map (Mask in Red)")
        plt.show()

    return M


def loadCESM1Data(fileNum: str, parameters: Dict[str, Any],
                  latIn=[0], lonIn=[0]) -> ClimateData:
    """Load CESM1 Dataset into Numpy arrays."""
    # pull variables from parameter dictionary
    runType = parameters['runType']
    baseDir = parameters['baseDir']
    dataSource = parameters['dataSource']
    dataFreq = parameters['dataFreq']
    CESM1variableName = parameters['CESM1variableName']
    ensembleNums = parameters['ensembleNums']
    start_dateInt = parameters['start_dateInt']
    stop_dateInt = parameters['stop_dateInt']
    isBulk = parameters['isBulk']
    hemisphere = parameters['hemisphere']

    if isBulk:
        CESM1variableName = CESM1variableName[4:]

    if runType == '40memberLargeEnsemble':

        # directory holding data for the 40-Member Large Ensemble
        # forcing data is contained in other data files,
        # so read forcings from ICEFRAC (T would also work)
        if CESM1variableName == 'forcing':
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                'ICEFRAC')
        else:
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                CESM1variableName)

        # sea ice thickness has separate data files
        # for the northern and southern hemispheres
        if CESM1variableName == 'HI':
            files = glob.glob(
                os.path.join(
                    dataDir,
                    '*.{:03d}.*.*{}h.*.nc'.format(
                        fileNum, hemisphere.lower())),
                recursive=True)

        elif CESM1variableName == 'flux':
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                'FSNT')
            files = glob.glob(
                os.path.join(
                    dataDir,
                    '*.{:03d}.*.nc'.format(fileNum)),
                recursive=True)

            dataDir2 = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                'FLNT')
            files2 = glob.glob(
                os.path.join(
                    dataDir2,
                    '*.{:03d}.*.nc'.format(fileNum)),
                recursive=True)
        
        else:
            glob_format_str = os.path.join(
                dataDir, '*.{:03d}.*.nc'.format(fileNum))
            files = glob.glob(glob_format_str, recursive=True)
        for ifile in range(len(files)):
            if CESM1variableName == 'ICEFRAC':
                # load in CESM1 ICEFRAC data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1ICEFRAC(files[ifile], 0)
            elif CESM1variableName == 'TS':
                # load in CESM1 TS data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1TS(files[ifile], 0)
            elif CESM1variableName == 'T':
                # load in CESM1 T data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1T(files[ifile], 0)

                # The observable T (atmospheric temperature)
                # contains data for multiple altitudes,
                # so select a single altitude for analysis
                # TODO: Revisit this!
                altInd = 0
                dataCur = np.squeeze(dataCur[:, :, :, altInd])

            elif CESM1variableName == 'SST':
                # load in CESM1 SST data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1SST(files[ifile], 0)

            elif CESM1variableName == 'TEMP':
                # load in CESM1 TEMP data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1TEMP(files[ifile], 0)

                # The observable TEMP (oceanic potential temperature)
                # contains data for multiple depths,
                # so select a single altitude for analysis
                # TODO: Revisit this!
                depthInd = 0
                dataCur = np.squeeze(dataCur[:, :, :, depthInd])

            elif CESM1variableName == 'HI':
                # load in CESM1 HI data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1HI(files[ifile], 0)

            elif CESM1variableName == 'flux':
                # load CESM1 FSNT and FLNT data,
                # and take their difference to get the net flux
                solar_flux, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1FSNT(files[ifile], 0)
                thermal_flux, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1FLNT(files2[ifile], 0)
                flux = solar_flux - thermal_flux
                dataCur = flux

            elif CESM1variableName == 'forcing':
                # load in CESM1 forcing data
                ch4vmr, co2vmr, f11vmr, f12vmr, n2ovmr, \
                    sol_tsi, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1forcing(files[ifile], 0)
                # Use only CO2 forcing for now.
                # We could also concatenate all of the forcings
                # into one data matrix, or load them all separately
                dataCur = co2vmr
                dataCur = np.expand_dims(dataCur, axis=(1, 2))

            # fix date offset
            dateIntCur[1:] = dateIntCur[:-1]
            dateIntCur[0] = dateIntCur[0]-100

            # initialize variables if first file being processed
            if ifile == 0:
                data = dataCur
                dateInt = dateIntCur
                lat = latCur
                lon = lonCur
                timeDays = timeDaysCur
            else:
                # other than first file, concatenate the data by time
                data = np.concatenate((data, dataCur), axis=0)
                dateInt = np.concatenate((dateInt, dateIntCur), axis=0)
                timeDays = np.concatenate((timeDays, timeDaysCur), axis=0)

        descStr = 'Ensemble Num {} {}'.format(fileNum, CESM1variableName)

    elif runType == 'AvgOfLargeEnsemble':
        # directory holding data for the Coupled Control Run
        # forcing data is contained in other data files,
        # so read forcings from ICEFRAC (T would also work)
        if CESM1variableName == 'forcing':
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                'ICEFRAC')
        else:
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                CESM1variableName)

        for iens in range(len(ensembleNums)):
            fileNum = ensembleNums[iens]

            # create list of all files from Ensemble Member fileNum
            # sea ice thickness has separate data files
            # for the northern and southern hemispheres
            if CESM1variableName == 'HI':
                files = glob.glob(
                    os.path.join(
                        dataDir,
                        '*.{:03d}.*.*{}h.*.nc'.format(
                            fileNum,
                            hemisphere.lower())),
                    recursive=True)

            elif CESM1variableName == 'flux':
                dataDir = os.path.join(
                    baseDir,
                    parameters['dataDirName'],
                    dataSource,
                    dataFreq,
                    'LargeEnsemble',
                    'FSNT')
                files = glob.glob(
                    os.path.join(
                        dataDir,
                        '*.{:03d}.*.nc'.format(fileNum)),
                    recursive=True)

                dataDir2 = os.path.join(
                    baseDir,
                    parameters['dataDirName'],
                    dataSource,
                    dataFreq,
                    'LargeEnsemble',
                    'FLNT')
                files2 = glob.glob(
                    os.path.join(
                        dataDir2,
                        '*.{:03d}.*.nc'.format(fileNum)),
                    recursive=True)

            else:
                files = glob.glob(
                    os.path.join(
                        dataDir,
                        '*.{:03d}.*.nc'.format(fileNum)),
                    recursive=True)

            for ifile in range(len(files)):
                if CESM1variableName == 'ICEFRAC':
                    # load in CESM1 ICEFRAC data
                    dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1ICEFRAC(files[ifile], 0)

                elif CESM1variableName == 'TS':
                    # load in CESM1 TS data
                    dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1TS(files[ifile], 0)

                elif CESM1variableName == 'T':
                    # load in CESM1 T data
                    dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1T(files[ifile], 0)

                    # The observable T (atmospheric temperature) contains data
                    # for multiple altitudes,
                    # so select a single altitude for analysis
                    # TODO: Revisit this!
                    altInd = 0
                    dataCur = np.squeeze(dataCur[:, :, :, altInd])

                elif CESM1variableName == 'SST':
                    # load in CESM1 SST data
                    dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1SST(files[ifile], 0)

                elif CESM1variableName == 'TEMP':
                    # load in CESM1 TEMP data
                    dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1TEMP(files[ifile], 0)

                    # The observable TEMP(oceanic potential temperature)
                    # contains data for multiple depths,
                    # so select a single altitude for analysis
                    # TODO: Revisit this!
                    depthInd = 0
                    dataCur = np.squeeze(dataCur[:, :, :, depthInd])
                elif CESM1variableName == 'HI':
                    # load in CESM1 HI data
                    dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1HI(files[ifile], 0)

                elif CESM1variableName == 'flux':
                    # load CESM1 FSNT and FLNT data,
                    # and take their difference to get the net flux
                    solar_flux, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1FSNT(files[ifile], 0)
                    thermal_flux, dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1FLNT(files2[ifile], 0)
                    flux = solar_flux - thermal_flux
                    dataCur = flux

                elif CESM1variableName == 'forcing':
                    # load in CESM1 forcing data
                    ch4vmr, co2vmr, f11vmr, f12vmr, n2ovmr, sol_tsi, \
                        dateIntCur, latCur, lonCur, timeDaysCur =\
                        readCESM1forcing(files[ifile], 0)

                    # Use only CO2 forcing for now.
                    # We could also concatenate all of the forcings
                    # into one data matrix, or load them all separately
                    dataCur = co2vmr
                    dataCur = np.expand_dims(dataCur, axis=(1, 2))

                # fix date offset
                dateIntCur[1:] = dateIntCur[:-1]
                dateIntCur[0] = dateIntCur[0]-100

                # initialize variables if first file being processed
                if ifile == 0:
                    data = dataCur
                    dateInt = dateIntCur
                    lat = latCur
                    lon = lonCur
                    timeDays = timeDaysCur
                else:
                    # other than first file, concatenate the data by time
                    data = np.concatenate((data, dataCur), axis=0)
                    dateInt = np.concatenate((dateInt, dateIntCur), axis=0)
                    timeDays = np.concatenate((timeDays, timeDaysCur), axis=0)

            # sort by date
            inds = np.argsort(dateInt)
            dateInt = np.sort(dateInt)
            data = data[inds, :, :]
            timeDays = timeDays[inds]

            # restrict data to time range covered
            # by all ensemble runs
            # #TODO: hardcoded
            t1 = 19200101
            t2 = 21001201
            ind1 = np.argwhere(dateInt == t1)[0, 0]
            ind2 = np.argwhere(dateInt == t2)[0, 0]
            data = data[ind1:ind2+1, :, :]

            # calculate cumulative moving average of data
            if iens == 0:
                data_all = data
                dateInt = dateInt[ind1:ind2+1]
                timeDays = timeDays[ind1:ind2+1]
            else:
                data_all = (data + iens*data_all)/(iens+1)

        data = data_all
        descStr = 'Ensemble Avg {}'.format(CESM1variableName)

    elif runType == 'CoupledControlRun':
        # directory holding data for the Coupled Control Run
        # forcing data is contained in other data files,
        # so read forcings from ICEFRAC (T would also work)
        if CESM1variableName == 'forcing':
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'CoupledControlRun',
                'ICEFRAC')
        else:
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'CoupledControlRun')

        # create list of all control run files
        # sea ice thickness has separate data files
        # for the northern and southern hemispheres
        if CESM1variableName == 'HI':
            files = glob.glob(
                os.path.join(
                    dataDir,
                    '*.*{}h.*.nc'.format(hemisphere.lower())),
                recursive=True)

        elif CESM1variableName == 'flux':
            dataDir = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                'FSNT')
            files = glob.glob(
                os.path.join(
                    dataDir,
                    '*.nc'),
                recursive=True)
            dataDir2 = os.path.join(
                baseDir,
                parameters['dataDirName'],
                dataSource,
                dataFreq,
                'LargeEnsemble',
                'FLNT')
            files2 = glob.glob(
                os.path.join(
                    dataDir2,
                    '*.nc'),
                recursive=True)
        else:
            files = glob.glob(
                os.path.join(dataDir, '*.nc'),
                recursive=True)
        files = glob.glob(
            os.path.join(dataDir, '*.nc'),
            recursive=True)

        for ifile in range(len(files)):
            if CESM1variableName == 'ICEFRAC':
                # load in CESM1 ICEFRAC data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1ICEFRAC(files[ifile], 0)

            elif CESM1variableName == 'TS':
                # load in CESM1 TS data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1TS(files[ifile], 0)

            elif CESM1variableName == 'T':
                # load in CESM1 T data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1T(files[ifile], 0)

                # The observable T (atmospheric temperature) contains data
                # for multiple altitudes,
                # so select a single altitude for analysis
                # TODO: Revisit this!
                altInd = 0
                dataCur = np.squeeze(dataCur[:, :, :, altInd])

            elif CESM1variableName == 'SST':
                # load in CESM1 SST data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1SST(files[ifile], 0)

            elif CESM1variableName == 'TEMP':
                # load in CESM1 TEMP data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1TEMP(files[ifile], 0)

                # The observable TEMP(oceanic potential temperature)
                # contains data for multiple depths,
                # so select a single altitude for analysis
                # TODO: Revisit this!
                depthInd = 0
                dataCur = np.squeeze(dataCur[:, :, :, depthInd])

            elif CESM1variableName == 'HI':
                # load in CESM1 HI data
                dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1HI(files[ifile], 0)

            elif CESM1variableName == 'flux':
                # load CESM1 FSNT and FLNT data,
                # and take their difference to get the net flux
                solar_flux, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1FSNT(files[ifile], 0)
                thermal_flux, dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1FLNT(files2[ifile], 0)
                flux = solar_flux - thermal_flux
                dataCur = flux

            elif CESM1variableName == 'forcing':
                # load in CESM1 forcing data
                ch4vmr, co2vmr, f11vmr, f12vmr, n2ovmr, sol_tsi, \
                    dateIntCur, latCur, lonCur, timeDaysCur =\
                    readCESM1forcing(files[ifile], 0)

                # Use only CO2 forcing for now.
                # We could also concatenate all of the forcings
                # into one data matrix, or load them all separately
                dataCur = co2vmr
                dataCur = np.expand_dims(dataCur, axis=(1, 2))

            # fix date offset
            dateIntCur[1:] = dateIntCur[:-1]
            dateIntCur[0] = dateIntCur[0]-100

            if ifile == 0:
                data = dataCur
                dateInt = dateIntCur
                lat = latCur
                lon = lonCur
                timeDays = timeDaysCur
            else:
                # other than first file, concatenate the data by time
                data = np.concatenate((data, dataCur), axis=0)
                dateInt = np.concatenate((dateInt, dateIntCur), axis=0)
                timeDays = np.concatenate((timeDays, timeDaysCur), axis=0)

        descStr = 'Coupled Control Run'

    # sort by date
    inds = np.argsort(dateInt)
    dateInt = np.sort(dateInt)
    data = data[inds, :, :]
    timeDays = timeDays[inds]

    # select time range of data
    timeInds = range(np.argwhere(dateInt == start_dateInt)[0, 0],
                     np.argwhere(dateInt == stop_dateInt)[0, 0] + 1)
    data = data[timeInds, :, :]
    dateInt = dateInt[timeInds]
    timeDays = timeDays[timeInds]

    # now lets interpolate any missing data
    # first we convert from 3d to 2d to match expected input
    # after the interpolation, we convert back
    lat_length = data.shape[2]
    lon_length = data.shape[1]
    data_2d = data.reshape(data.shape[0], -1)
    interpolator = Interpolator()
    data_2d, dateInt = interpolator.interpolate_missing_data(data_2d, dateInt)
    data = data_2d.reshape(
        data_2d.shape[0],
        lon_length,
        lat_length
    )

    # round sea ice concentration value to the nearest percent
    if CESM1variableName == 'ICEFRAC':
        data[data > 100.0] = 100.0
        data[data < 0.0] = 0.0

    # cutoff temperature placeholders
    if (CESM1variableName == 'TEMP' or
            CESM1variableName == 'SST' or
            CESM1variableName == 'T' or
            CESM1variableName == 'TS' or
            CESM1variableName == 'HI'):
        data[data > 1e20] = 0.0

    # round to nearest percent
    if parameters['roundConc']:
        data = np.around(data)

    # sea ice concentration percentages
    # at or below this value will be set to zero
    if (CESM1variableName == 'ICEFRAC' and
            parameters['dataSource'] != 'err'):
        data[data <= parameters['concentrationToConsiderZero']] = 0.0

    # The TEMP,
    # SST (ocean potential temperature),
    # and HI (sea ice thickness)
    # variables are given on a 2D lat/lon grid,
    # so interpolate using input lat/lon array values
    # to make it consistent with other variables like ICEFRAC
    if (CESM1variableName == 'TEMP' or
            CESM1variableName == 'SST' or
            CESM1variableName == 'HI'):
        xx, yy = np.meshgrid(lonIn, latIn)
        lat1D = np.reshape(lat, (lat.size))
        lon1D = np.reshape(lon, (lon.size))
        dataInterp = np.zeros((len(data), len(latIn), len(lonIn)), dtype=float)

        # TODO: CURRENTLY A HUGE PERFORMANCE BOTTLENECK
        # loop over first dimension of data, which is time
        for itime in range(len(data)):
            dataCur = np.squeeze(data[itime, :, :])
            dataCur1D = np.reshape(dataCur, (dataCur.size))
            dataInterp[itime, :, :] = interpolate.griddata(
                (lon1D, lat1D), dataCur1D, (xx, yy), method='nearest')

        data = dataInterp
        lat = latIn
        lon = lonIn

    if isBulk:
        descStr = descStr + ' (bulk)'

    return ClimateData(data, dateInt, lat, lon, timeDays, descStr)


def loadNSIDCData(parameters: Dict[str, Any],
                  latIn=[0], lonIn=[0]) -> ClimateData:
    """Load NSIDC Dataset into Numpy arrays."""
    # pull variables from parameter dictionary
    baseDir = parameters['baseDir']
    dataSource = parameters['dataSource']
    dataFreq = parameters['dataFreq']
    NSIDCvariableName = parameters['NSIDCvariableName']
    hemisphere = parameters['hemisphere']
    start_dateInt = parameters['start_dateInt']
    stop_dateInt = parameters['stop_dateInt']
    flag_use_missing = parameters['flag_use_missing']
    isBulk = parameters['isBulk']

    if isBulk:
        NSIDCvariableName = NSIDCvariableName[4:]

    if NSIDCvariableName == 'merged':
        if hemisphere == 'N':
            dataSetStr = 'north_v4_regrid'
        elif hemisphere == 'S':
            dataSetStr = 'south_v4_regrid'

    elif NSIDCvariableName == 'sst':
        dataSetStr = 'ORAS5'

    elif NSIDCvariableName == 't2m':
        dataSetStr = 'ERA5'

    elif NSIDCvariableName == 'flux':
        dataSetStr = 'ERA5'

    elif NSIDCvariableName == 'thickness':
        dataSetStr = 'PIOMAS'

    # directory holding data for the 40-Member Large Ensemble
    dataDir = os.path.join(
        baseDir,
        parameters['dataDirName'],
        dataSource,
        dataFreq,
        dataSetStr,
        NSIDCvariableName)

    # create list of all observational data files.
    # Note that flux observable requires loading
    # two different observables (tsr and ttr)
    files = glob.glob(
        os.path.join(dataDir, '*.nc'),
        recursive=True)

    for ifile in range(len(files)):
        if NSIDCvariableName == 'merged':
            # load NSIDC data
            dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                readNSIDCSIC(files[ifile], NSIDCvariableName, 0)

        elif NSIDCvariableName == 'sst':
            # load ORAS5 data
            dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                read_ORAS5_sst(files[ifile], 0)

        elif NSIDCvariableName == 't2m':
            # load ERA5 data
            dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                read_ERA5_t2m(files[ifile], 0)

        elif NSIDCvariableName == 'top_flux':
            # load ERA5 top flux data
            dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
               read_ERA5_top_flux(files[ifile], 0)

        elif NSIDCvariableName == 'surf_flux':
            # load ERA5 surface flux data
            dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
               read_ERA5_surf_flux(files[ifile], 0)

        elif NSIDCvariableName == 'thickness':
            # load thickness data
            dataCur, dateIntCur, latCur, lonCur, timeDaysCur =\
                read_PSC_thickness(files[ifile], 0)

        if ifile == 0:
            data = dataCur
            dateInt = dateIntCur
            lat = latCur
            lon = lonCur
            timeDays = timeDaysCur
        else:
            # other than first file, concatenate the data by time
            data = np.concatenate((data, dataCur), axis=0)
            dateInt = np.concatenate((dateInt, dateIntCur), axis=0)
            timeDays = np.concatenate((timeDays, timeDaysCur), axis=0)
    
    # sort by date
    inds = np.argsort(dateInt)
    dateInt = np.sort(dateInt)
    data = data[inds, :, :]
    timeDays = timeDays[inds]

    # select time range of data
    timeInds = range(np.argwhere(dateInt == start_dateInt)[0, 0],
                     np.argwhere(dateInt == stop_dateInt)[0, 0] + 1)
    data = data[timeInds, :, :]
    dateInt = dateInt[timeInds]
    timeDays = timeDays[timeInds]
    
    # now lets interpolate any missing data
    # first we convert from 3d to 2d to match expected input
    # after the interpolation, we convert back
    lat_length = data.shape[2]
    lon_length = data.shape[1]
    data_2d = data.reshape(data.shape[0], -1)
    interpolator = Interpolator()
    data_2d, dateInt = interpolator.interpolate_missing_data(data_2d, dateInt)
    data = data_2d.reshape(
        data_2d.shape[0],
        lon_length,
        lat_length
    )

    descStr = 'NSIDC {}'.format(NSIDCvariableName)

    # round sea ice concentration value to the nearest percent
    if NSIDCvariableName == 'merged':
        data[data > 100.0] = 100.0
        data[data < 0.0] = 0.0

    # cutoff temperature placeholders
    if NSIDCvariableName == 'sst' or NSIDCvariableName == 't2m':
        data[data > 1e20] = 0.0

    # round to nearest percent
    if parameters['roundConc']:
        data = np.around(data)

    # sea ice concentration percentages at or
    # below this value will be set to zero
    if (NSIDCvariableName == 'merged' and
            parameters['dataSource'] != 'err'):
        data[data <= parameters['concentrationToConsiderZero']] = 0

    # The ORAS5, ERA5, and PSC PIOMAS thickness data
    # have different lat and lon arrays than the NSIDC data,
    # so interpolate using the NSIDC lat/lon array values
    # to make them consistent with the NSIDC data
    if parameters['NSIDCvariableName'][-len('merged'):] != 'merged':
        if NSIDCvariableName == 'thickness':
            xx, yy = np.meshgrid(lonIn, latIn)
            lat1D = np.reshape(lat, (lat.size))
            lon1D = np.reshape(lon, (lon.size))
            dataInterp = np.zeros((len(data), len(latIn), len(lonIn)),
                                  dtype=float)

            # # TODO: CURRENTLY A HUGE PERFORMANCE BOTTLENECK
            # loop over first dimension of data, which is time
            for itime in range(len(data)):
                dataCur = np.squeeze(data[itime, :, :])
                dataCur1D = np.reshape(dataCur, (dataCur.size))
                dataInterp[itime, :, :] = interpolate.griddata(
                    (lon1D, lat1D), dataCur1D, (xx, yy),
                    method='nearest')

            data = dataInterp
            lat = latIn
            lon = lonIn

        else:
            dataInterp = np.zeros((len(data), len(latIn), len(lonIn)),
                                  dtype=float)

            # loop over first dimension of data, which is time
            for itime in range(len(data)):
                f = interpolate.interp2d(
                    lon, lat, np.squeeze(data[itime, :, :]))
                dataInterp[itime, :, :] = f(lonIn, latIn)
            data = dataInterp
            lat = latIn
            lon = lonIn

    if isBulk:
        descStr = descStr + ' (bulk)'

    return ClimateData(data, dateInt, lat, lon, timeDays, descStr)


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
    temp_parameters = {}
    temp_parameters['baseDir'] = parameters['baseDir']
    temp_parameters['dataDirName'] = parameters['dataDirName']
    temp_parameters['dataSource'] = parameters['dataSource']
    temp_parameters['dataFreq'] = parameters['dataFreq']
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
    climate_data = loadNSIDCData(temp_parameters)
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
