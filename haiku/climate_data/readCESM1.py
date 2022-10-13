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

from typing import Tuple
import numpy as np
import numpy.typing as npt
from netCDF4 import Dataset
import calendar
import datetime


def readCESM1ICEFRAC(fnFull: str,
                     VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                             npt.NDArray, npt.NDArray,
                                             npt.NDArray]:
    """Read ICEFRAC (sea ice concentration) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # ICEFRAC
    #            Size:       Tx192x288
    #            Dimensions: time,lat,lon
    #            Datatype:   single
    #            Attributes: units        = 'fraction'
    #            long_name    = 'Fraction of sfc area covered by sea-ice'
    # Convert ICEFRAC from fraction to percent
    icefrac = 100.0 * np.array(ncfile.variables['ICEFRAC'])

    # date
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   int32
    #            Attributes: long_name = 'current date (YYYYMMDD)'
    dateInt = np.array(ncfile.variables['date'])

    # lat
    #            Size:       192x1
    #            Dimensions: lat
    #            Datatype:   double
    #            Attributes: long_name = 'latitude'
    #                        units     = 'degrees_north'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #            Size:       288x1
    #            Dimensions: lon
    #            Datatype:   double
    #            Attributes: long_name = 'longitude'
    #                        units     = 'degrees_east'
    lon = np.array(ncfile.variables['lon'])

    # time
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   double
    #            Attributes: long_name = 'time'
    #                        units     = 'days since 1920-01-01 00:00:00'
    #                        calendar  = 'noleap'
    timeDays = np.array(ncfile.variables['time'])

    ncfile.close()

    return icefrac, dateInt, lat, lon, timeDays


def readCESM1TS(fnFull: str,
                VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                        npt.NDArray, npt.NDArray,
                                        npt.NDArray]:
    """Read TS (atmospheric surface temperature) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # TS
    #        Size:       Tx192x288
    #        Dimensions: time,lat,lon
    #        Datatype:   single
    #        Attributes:
    #                    mdims        = 1
    #                    units        = 'K'
    #                    long_name    = 'Surface temperature (radiative)'
    #                    cell_methods = 'time: mean'
    TS = np.array(ncfile.variables['TS'])

    # date
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   int32
    #            Attributes: long_name = 'current date (YYYYMMDD)'
    dateInt = np.array(ncfile.variables['date'])

    # lat
    #            Size:       192x1
    #            Dimensions: lat
    #            Datatype:   double
    #            Attributes: long_name = 'latitude'
    #                        units     = 'degrees_north'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #            Size:       288x1
    #            Dimensions: lon
    #            Datatype:   double
    #            Attributes: long_name = 'longitude'
    #                        units     = 'degrees_east'
    lon = np.array(ncfile.variables['lon'])

    # time
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   double
    #            Attributes: long_name = 'time'
    #                        units     = 'days since 1920-01-01 00:00:00'
    #                        calendar  = 'noleap'
    timeDays = np.array(ncfile.variables['time'])

    ncfile.close()

    return TS, dateInt, lat, lon, timeDays


def readCESM1T(fnFull: str,
               VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                       npt.NDArray, npt.NDArray,
                                       npt.NDArray]:
    """Read T (atmospheric temperature) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # T
    #        Size:       Tx30x192x288
    #        Dimensions: time,lev,lat,lon
    #        Datatype:   single
    #        Attributes:
    #                    mdims        = 1
    #                    units        = 'K'
    #                    long_name    = 'Temperature'
    #                    cell_methods = 'time: mean'
    # Permute dimensions to match the dimension order in ICEFRAC
    T = np.transpose(np.array(ncfile.variables['T']), (0, 2, 3, 1))

    # date
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   int32
    #            Attributes: long_name = 'current date (YYYYMMDD)'
    dateInt = np.array(ncfile.variables['date'])

    # lat
    #            Size:       192x1
    #            Dimensions: lat
    #            Datatype:   double
    #            Attributes: long_name = 'latitude'
    #                        units     = 'degrees_north'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #            Size:       288x1
    #            Dimensions: lon
    #            Datatype:   double
    #            Attributes: long_name = 'longitude'
    #                        units     = 'degrees_east'
    lon = np.array(ncfile.variables['lon'])

    # time
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   double
    #            Attributes: long_name = 'time'
    #                        units     = 'days since 1920-01-01 00:00:00'
    #                        calendar  = 'noleap'
    timeDays = np.array(ncfile.variables['time'])

    ncfile.close()

    return T, dateInt, lat, lon, timeDays


def readCESM1SST(fnFull: str,
                 VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                         npt.NDArray, npt.NDArray,
                                         npt.NDArray]:
    """Read SST (sea surface temperature) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # SST
    #        Size:       Tx1x384x320
    #        Dimensions: time,z_t,nlat,nlon
    #        Datatype:   single
    #        Attributes:
    #                    long_name     = 'Potential Temperature'
    #                    units         = 'degC'
    #                    coordinates   = 'TLONG TLAT z_t time'
    #                    grid_loc      = '3111'
    #                    cell_methods  = 'time: mean'
    #                    _FillValue    = 9.969209968386869e+36
    #                    missing_value = 9.969209968386869e+36
    SST = np.array(ncfile.variables['SST'])
    # Permute dimensions to match the dimension order in ICEFRAC
    SST = np.transpose(SST, (0, 2, 3, 1))

    # TLAT
    #        Size:       320x384
    #        Dimensions: nlon,nlat
    #        Datatype:   double
    #        Attributes:
    #                    long_name     = 'array of t-grid latitudes'
    #                    units         = 'degrees_north'
    #                    _FillValue    = 9.969209968386869e+36
    #                    missing_value = 9.969209968386869e+36
    lat = np.array(ncfile.variables['TLAT'])

    # TLONG
    #        Size:       320x384
    #        Dimensions: nlon,nlat
    #        Datatype:   double
    #        Attributes:
    #                    long_name     = 'array of t-grid longitudes'
    #                    units         = 'degrees_east'
    #                    _FillValue    = 9.969209968386869e+36
    #                    missing_value = 9.969209968386869e+36
    lon = np.array(ncfile.variables['TLONG'])

    # time
    #        Size:       Tx1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'time'
    #                    units     = 'days since 0000-01-01 00:00:00'
    #                    bounds    = 'time_bound'
    #                    calendar  = 'noleap'
    timeDays = np.array(ncfile.variables['time']).astype(int)

    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])+1
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0]) + extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1
        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m%d'))

    ncfile.close()

    return SST, dateInt, lat, lon, timeDays


def readCESM1TEMP(fnFull: str,
                  VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                          npt.NDArray, npt.NDArray,
                                          npt.NDArray]:
    """Read TEMP (oceanic potential temperature) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # TEMP
    #        Size:       Tx60x384x320
    #        Dimensions: time,z_t,nlat,nlon
    #        Datatype:   single
    #        Attributes:
    #                    long_name     = 'Potential Temperature'
    #                    units         = 'degC'
    #                    coordinates   = 'TLONG TLAT z_t time'
    #                    grid_loc      = '3111'
    #                    cell_methods  = 'time: mean'
    #                    _FillValue    = 9.969209968386869e+36
    #                    missing_value = 9.969209968386869e+36
    TEMP = np.array(ncfile.variables['TEMP'])
    # Permute dimensions to match the dimension order in ICEFRAC
    TEMP = np.transpose(TEMP, (0, 2, 3, 1))

    # TLAT
    #        Size:       320x384
    #        Dimensions: nlon,nlat
    #        Datatype:   double
    #        Attributes:
    #                    long_name     = 'array of t-grid latitudes'
    #                    units         = 'degrees_north'
    #                    _FillValue    = 9.969209968386869e+36
    #                    missing_value = 9.969209968386869e+36
    lat = np.array(ncfile.variables['TLAT'])

    # TLONG
    #        Size:       320x384
    #        Dimensions: nlon,nlat
    #        Datatype:   double
    #        Attributes:
    #                    long_name     = 'array of t-grid longitudes'
    #                    units         = 'degrees_east'
    #                    _FillValue    = 9.969209968386869e+36
    #                    missing_value = 9.969209968386869e+36
    lon = np.array(ncfile.variables['TLONG'])

    # time
    #        Size:       Tx1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'time'
    #                    units     = 'days since 0000-01-01 00:00:00'
    #                    bounds    = 'time_bound'
    #                    calendar  = 'noleap'
    timeDays = np.array(ncfile.variables['time']).astype(int)

    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])+1
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0]) + extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1
        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m%d'))

    ncfile.close()

    return TEMP, dateInt, lat, lon, timeDays


def readCESM1HI(fnFull: str,
                VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                        npt.NDArray, npt.NDArray,
                                        npt.NDArray]:
    """Read HI (sea ice thickness) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # hi
    #        Size:       320x{104,76}x1872
    #        Dimensions: ni,nj,time
    #        Datatype:   single
    #        Attributes:
    #                    units         = 'm'
    #                    long_name     = 'grid cell mean ice thickness'
    #                    coordinates   = 'TLON TLAT time'
    #                    cell_measures = 'area: tarea'
    #                    missing_value = 1.000000015047466e+30
    #                    _FillValue    = 1.000000015047466e+30
    #                    comment       = 'ice volume per unit grid cell area'
    #                    cell_methods  = 'time: mean'
    #                    time_rep      = 'averaged'
    THICKNESS = np.array(ncfile.variables['hi'])

    # TLAT
    #        Size:       320x{104,76}
    #        Dimensions: ni,nj
    #        Datatype:   single
    #        Attributes:
    #                    long_name     = 'T grid center latitude'
    #                    units         = 'degrees_north'
    #                    missing_value = 1.000000015047466e+30
    #                    _FillValue    = 1.000000015047466e+30
    #                    bounds        = 'latt_bounds'
    lat = np.array(ncfile.variables['TLAT'])

    # TLON
    #        Size:       320x{104,76}
    #        Dimensions: ni,nj
    #        Datatype:   single
    #        Attributes:
    #                    long_name     = 'T grid center longitude'
    #                    units         = 'degrees_east'
    #                    missing_value = 1.000000015047466e+30
    #                    _FillValue    = 1.000000015047466e+30
    #                    bounds        = 'lont_bounds'
    lon = np.array(ncfile.variables['TLON'])

    # time
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   single
    #        Attributes:
    #                    long_name = 'model time'
    #                    units     = 'days since 0000-01-01 00:00:00'
    #                    calendar  = 'noleap'
    #                    bounds    = 'time_bounds'
    timeDays = np.array(ncfile.variables['time']).astype(int)

    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])+1
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0]) + extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1
        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m%d'))

    ncfile.close()

    return THICKNESS, dateInt, lat, lon, timeDays


def readCESM1FSNT(fnFull: str,
                  VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                          npt.NDArray, npt.NDArray,
                                          npt.NDArray]:
    """Read FSNT (net solar flux at top of model) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # FSNT
    #        Size:       288x192x1872
    #        Dimensions: lon,lat,time
    #        Datatype:   single
    #        Attributes:
    #                    Sampling_Sequence = 'rad_lwsw'
    #                    units             = 'W/m2'
    #                    long_name         = 'Net solar flux at top of model'
    #                    cell_methods      = 'time: mean'
    solar_flux = np.array(ncfile.variables['FSNT'])

    # lat
    #        Size:       192x1
    #        Dimensions: lat
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'latitude'
    #                    units     = 'degrees_north'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #        Size:       288x1
    #        Dimensions: lon
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'longitude'
    #                    units     = 'degrees_east'
    lon = np.array(ncfile.variables['lon'])

    # time
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'time'
    #                    units     = 'days since 1850-01-01 00:00:00'
    #                    calendar  = 'noleap'
    #                    bounds    = 'time_bnds'
    timeDays = np.array(ncfile.variables['time']).astype(int)

    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])+1
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0]) + extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1
        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m%d'))

    ncfile.close()

    return solar_flux, dateInt, lat, lon, timeDays


def readCESM1FLNT(fnFull: str,
                  VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                          npt.NDArray, npt.NDArray,
                                          npt.NDArray]:
    """Read FLNT (net longwave flux at top of model) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # FLNT
    #        Size:       288x192x1872
    #        Dimensions: lon,lat,time
    #        Datatype:   single
    #        Attributes:
    #           Sampling_Sequence = 'rad_lwsw'
    #           units             = 'W/m2'
    #           long_name         = 'Net longwave flux at top of model'
    #           cell_methods      = 'time: mean'
    thermal_flux = np.array(ncfile.variables['FLNT'])

    # lat
    #        Size:       192x1
    #        Dimensions: lat
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'latitude'
    #                    units     = 'degrees_north'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #        Size:       288x1
    #        Dimensions: lon
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'longitude'
    #                    units     = 'degrees_east'
    lon = np.array(ncfile.variables['lon'])

    # time
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'time'
    #                    units     = 'days since 1850-01-01 00:00:00'
    #                    calendar  = 'noleap'
    #                    bounds    = 'time_bnds'
    timeDays = np.array(ncfile.variables['time']).astype(int)

    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])+1
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0]) + extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1
        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m%d'))

    ncfile.close()

    return thermal_flux, dateInt, lat, lon, timeDays


def readCESM1FSURF_AI(fnFull: str,
                      VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                              npt.NDArray, npt.NDArray,
                                              npt.NDArray]:
    """Read FSURF_AI (net surface heat flux) values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # FSURF_AI
    #        Size:       288x192x1872
    #        Dimensions: lon,lat,time
    #        Datatype:   single
    #        Attributes:
    #           Sampling_Sequence = 'rad_lwsw'
    #           units             = 'W/m2'
    #           long_name         = 'Net longwave flux at top of model'
    #           cell_methods      = 'time: mean'
    flux = np.array(ncfile.variables['FSURF_AI'])

    # lat
    #        Size:       192x1
    #        Dimensions: lat
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'latitude'
    #                    units     = 'degrees_north'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #        Size:       288x1
    #        Dimensions: lon
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'longitude'
    #                    units     = 'degrees_east'
    lon = np.array(ncfile.variables['lon'])

    # time
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'time'
    #                    units     = 'days since 1850-01-01 00:00:00'
    #                    calendar  = 'noleap'
    #                    bounds    = 'time_bnds'
    timeDays = np.array(ncfile.variables['time']).astype(int)

    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])+1
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0]) + extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1
        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m%d'))

    ncfile.close()

    return flux, dateInt, lat, lon, timeDays


def readCESM1forcing(fnFull: str,
                     VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                             npt.NDArray, npt.NDArray,
                                             npt.NDArray, npt.NDArray,
                                             npt.NDArray, npt.NDArray,
                                             npt.NDArray, npt.NDArray]:
    """Read greenhouse gas forcing values.

    Input: CESM1 LENS data files in NetCDF format.
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # ch4vmr
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'ch4 volume mixing ratio'
    ch4vmr = np.array(ncfile.variables['ch4vmr'])

    # co2vmr
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'co2 volume mixing ratio'
    co2vmr = np.array(ncfile.variables['co2vmr'])

    # f11vmr
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'f11 volume mixing ratio'
    f11vmr = np.array(ncfile.variables['f11vmr'])

    # f12vmr
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'f12 volume mixing ratio'
    f12vmr = np.array(ncfile.variables['f12vmr'])

    # n2ovmr
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'n2o volume mixing ratio'
    n2ovmr = np.array(ncfile.variables['n2ovmr'])

    # sol_tsi
    #        Size:       1872x1
    #        Dimensions: time
    #        Datatype:   double
    #        Attributes:
    #                    long_name = 'total solar irradiance'
    #                    units     = 'W/m2'
    sol_tsi = np.array(ncfile.variables['sol_tsi'])

    # date
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   int32
    #            Attributes: long_name = 'current date (YYYYMMDD)'
    dateInt = np.array(ncfile.variables['date'])

    # lat
    #            Size:       192x1
    #            Dimensions: lat
    #            Datatype:   double
    #            Attributes: long_name = 'latitude'
    #                        units     = 'degrees_north'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #            Size:       288x1
    #            Dimensions: lon
    #            Datatype:   double
    #            Attributes: long_name = 'longitude'
    #                        units     = 'degrees_east'
    lon = np.array(ncfile.variables['lon'])

    # time
    #            Size:       Tx1
    #            Dimensions: time
    #            Datatype:   double
    #            Attributes: long_name = 'time'
    #                        units     = 'days since 1920-01-01 00:00:00'
    #                        calendar  = 'noleap'
    timeDays = np.array(ncfile.variables['time'])

    ncfile.close()

    return ch4vmr, co2vmr, f11vmr, f12vmr, \
        n2ovmr, sol_tsi, dateInt, lat, lon, timeDays
