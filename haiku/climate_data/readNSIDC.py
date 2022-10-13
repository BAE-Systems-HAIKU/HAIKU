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
import os


def readNSIDCSIC(fnFull: str, NSIDCvariableName: str,
                 VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                         npt.NDArray, npt.NDArray,
                                         npt.NDArray]:
    """Read sea ice concentration values from NSIDC files in NetCDF format.

    Importantly, the final sea ice concentration values are created
    by combining measurements from different sensors.
    There are different methods to produce slightly different final results:
        nsidc replaces goddard in V4
        goddard_merged_seaice_conc_monthly
        goddard_nt_seaice_conc_monthly
        goddard_bt_seaice_conc_monthly
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # goddard_merged_seaice_conc_monthly
    #            Size:       1x192x288
    #            Dimensions: time,lat,lon
    #            Datatype:   int16
    #            Attributes: standard_name: sea_ice_area_fraction
    #                        long_name: Goddard Edited Climate Data Record of
    #                        Passive Microwave Monthly Northern Hemisphere Sea
    #                        Ice Concentration, Goddard Edited
    #                        units: 1
    #                        add_offset: 0.0
    #                        scale_factor: 0.01
    #                        _FillValue: -1
    #                        missing_value: -1
    #                        flag_values: [-5 -4 -3 -2 -1]
    #                        flag_meanings: pole_hole lakes
    #                                       coastal land_mask missing_data
    #                        datum: +ellps=urn:ogc:def:crs:EPSG::4326
    #                        reference: http://nsidc.org/data/nsidc-0051.html
    #                        http://nsidc.org/data/nsidc-0079.html

    # Convert ICEFRAC from fraction to percent
    SIC_merged = 100.0*np.array(
        #v3 only; ncfile.variables['goddard_merged_seaice_conc_monthly'])
        ncfile.variables['cdr_seaice_conc_monthly'])

    # goddard_nt_seaice_conc_monthly
    #            Size:       1x192x288
    #            Dimensions: time,lat,lon
    #            Datatype:   int16
    #            Attributes: standard_name: sea_ice_area_fraction
    #                        long_name: Passive Microwave Monthly Northern
    #                                   Hemisphere Sea Ice Concentration by
    #                                   NASA Team algorithm with Goddard QC
    #                        units: 1
    #                        add_offset: 0.0
    #                        scale_factor: 0.01
    #                        _FillValue: -1
    #                        missing_value: -1
    #                        flag_values: [-5 -4 -3 -2 -1]
    #                        flag_meanings: pole_hole unused coastal
    #                                       land_mask missing_data
    #                        datum: +ellps=urn:ogc:def:crs:EPSG::4326
    #                        reference: Documentation available at:
    #                        http://nsidc.org/data/nsidc-0051.html

    # Convert ICEFRAC from fraction to percent
    SIC_nt = 100.0*np.array(ncfile.variables['nsidc_nt_seaice_conc_monthly'])

    # goddard_bt_seaice_conc_monthly
    #            Size:       1x192x288
    #            Dimensions: time,lat,lon
    #            Datatype:   int16
    #            Attributes: standard_name: sea_ice_area_fraction
    #                        long_name: Passive Microwave Monthly Northern
    #                                   Hemisphere Sea Ice Concentration by
    #                                   Bootstrap algorithm with Goddard QC
    #                        units: 1
    #                        add_offset: 0.0
    #                        scale_factor: 0.01
    #                        _FillValue: -1
    #                        missing_value: -1
    #                        flag_values: [-5 -4 -3 -2 -1]
    #                        flag_meanings: pole_hole unused coastal
    #                                       land_mask missing_data
    #                        datum: +ellps=urn:ogc:def:crs:EPSG::4326
    #                        reference: Documentation available at:
    #                                   http://nsidc.org/data/nsidc-0079.html

    # Convert ICEFRAC from fraction to percent
    SIC_bt = 100.0*np.array(ncfile.variables['nsidc_bt_seaice_conc_monthly'])

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
    #            Size:       1x1
    #            Dimensions: time
    #            Datatype:   double
    #            Attributes: long_name: ANSI date
    #                        units: days since 1601-01-01 00:00:00
    #                        calendar: gregorian
    timeDays = np.array(ncfile.variables['time'])

    ncfile.close()

    if NSIDCvariableName == 'merged':
        SIC = SIC_merged
    elif NSIDCvariableName == 'nt':
        SIC = SIC_nt
    elif NSIDCvariableName == 'bt':
        SIC = SIC_bt

    # instead of grabbing the date from the filename,
    # we can compute it from time variable

    # Let's calculate the dateInt string from the filename,
    # to support the .nc files we created to fill
    # in the missing months of NSIDC data
    #TODO: extract data from .nc file, not from filename
    dateInt = [int(fnFull.split('_')[-3] + '01')]

    # replace nan values with 100.0 in case of any MATLAB generated files
    SIC[np.isnan(SIC)] = 0.0

    # print( np.amax(SIC) )

    return SIC, dateInt, lat, lon, timeDays


def read_ERA5_t2m(fnFull: str,
                  VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                          npt.NDArray, npt.NDArray,
                                          npt.NDArray]:
    """Read 2m temperature values from ERA5 atmospheric forcing reanalysis.

    https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)
    # t2m
    #        Size:       132x721x1440
    #        Dimensions: time, latitude, longitude
    #        Datatype: ncfile.variables  int16
    #        Attributes:
    #                    scale_factor  = 0.0013174
    #                    add_offset    = 266.7527
    #                    _FillValue    = -32767
    #                    missing_value = -32767
    #                    units         = 'K'
    #                    long_name     = '2 metre temperature'
    t2m = np.array(ncfile.variables['t2m'])

    # latitude
    #        Size:       721x1
    #        Dimensions: latitude
    #        Datatype:   single
    #        Attributes:
    #                    units     = 'degrees_north'
    #                    long_name = 'latitude'
    lat = np.array(ncfile.variables['latitude'])

    # longitude
    #        Size:       1440x1
    #        Dimensions: longitude
    #        Datatype:   single
    #        Attributes:
    #                    units     = 'degrees_east'
    #                    long_name = 'longitude'
    lon = np.array(ncfile.variables['longitude'])

    # time
    #        Size:       1x1
    #        Dimensions: time
    #        Datatype:   int32
    #        Attributes:
    #                    units     = 'hours since 1900-01-01 00:00:00.0'
    #                    long_name = 'time'
    #                    calendar  = 'gregorian'
    timeDays = np.array(ncfile.variables['time'])/24
    
    ncfile.close()

    # Let's calculate the dateInt string from the filename
    # (note that these filenames were set by hand after
    # downloading from the Copernicus Climate Data Store)
    dateInt = np.empty_like(timeDays,dtype=int)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0]) + extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1

        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m') + '01')
    return t2m, dateInt, lat, lon, timeDays


def read_ERA5_top_flux(fnFull: str,
                       VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                               npt.NDArray, npt.NDArray,
                                               npt.NDArray]:
    """Read top net solar radiation and top net thermal radiation.

    Input: From the ERA5 atmospheric forcing reanalysis data in NetCDF format.
    Output: Net flux.

    https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # tsr
    #        Size:       1440x721xnTimes
    #        Dimensions: longitude,latitude,time
    #        Datatype:   int16
    #        Attributes:
    #                    scale_factor  = 606.3715
    #                    add_offset    = 19868368.8142
    #                    _FillValue    = -32767
    #                    missing_value = -32767
    #                    units         = 'J m**-2'
    #                    long_name     = 'Top net solar radiation'
    #                    standard_name = 'toa_net_upward_shortwave_flux'
    tsr = np.array(ncfile.variables['tsr'])

    # ttr
    #        Size:       1440x721xnTimes
    #        Dimensions: longitude,latitude,time
    #        Datatype:   int16
    #        Attributes:
    #                    scale_factor  = 341.5918
    #                    add_offset    = -19230458.7959
    #                    _FillValue    = -32767
    #                    missing_value = -32767
    #                    units         = 'J m**-2'
    #                    long_name     = 'Top net thermal radiation'
    #                    standard_name = 'toa_outgoing_longwave_flux'
    ttr = np.array(ncfile.variables['ttr'])

    # latitude
    #        Size:       721x1
    #        Dimensions: latitude
    #        Datatype:   single
    #        Attributes:
    #                    units     = 'degrees_north'
    #                    long_name = 'latitude'
    lat = np.array(ncfile.variables['latitude'])

    # longitude
    #        Size:       1440x1
    #        Dimensions: longitude
    #        Datatype:   single
    #        Attributes:
    #                    units     = 'degrees_east'
    #                    long_name = 'longitude'
    lon = np.array(ncfile.variables['longitude'])

    # time
    #        Size:       nTimes x 1
    #        Dimensions: time
    #        Datatype:   int32
    #        Attributes:
    #                    units     = 'hours since 1900-01-01 00:00:00.0'
    #                    long_name = 'time'
    #                    calendar  = 'gregorian'
    timeDays = np.array(ncfile.variables['time'])/24

    ncfile.close()

    # Let's calculate the dateInt string from the filename
    # (note that these filenames were
    # set by hand after downloading from the Copernicus Climate Data Store)
    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0])+extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1

        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m') + '01')

    # net flux is sum of tsr and ttr (due to sign conventions),
    # and divide by seconds/day to convert from J/m^2 to W/m^2
    # because tsr and ttr are accumulations scaled to be over one day
    flux = (tsr + ttr) / 86400

    return flux, dateInt, lat, lon, timeDays


def read_ERA5_surf_flux(fnFull: str,
                        VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                                npt.NDArray, npt.NDArray,
                                                npt.NDArray]:
    """Read surface net solar radiation and surface net thermal radiation.

    Input: ERA5 atmospheric forcing reanalysis data in NetCDF format.
    Output: Net flux.

    https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # ssr
    #        Size:       1440x721xnTimes
    #        Dimensions: longitude,latitude,time
    #        Datatype:   int16
    #        Attributes:
    #                    scale_factor  = 606.3715
    #                    add_offset    = 19868368.8142
    #                    _FillValue    = -32767
    #                    missing_value = -32767
    #                    units         = 'J m**-2'
    #                    long_name     = 'Top net solar radiation'
    #                    standard_name = 'toa_net_upward_shortwave_flux'
    ssr = np.array(ncfile.variables['ssr'])

    # str
    #        Size:       1440x721xnTimes
    #        Dimensions: longitude,latitude,time
    #        Datatype:   int16
    #        Attributes:
    #                    scale_factor  = 341.5918
    #                    add_offset    = -19230458.7959
    #                    _FillValue    = -32767
    #                    missing_value = -32767
    #                    units         = 'J m**-2'
    #                    long_name     = 'Top net thermal radiation'
    #                    standard_name = 'toa_outgoing_longwave_flux'
    str = np.array(ncfile.variables['str'])

    # latitude
    #        Size:       721x1
    #        Dimensions: latitude
    #        Datatype:   single
    #        Attributes:
    #                    units     = 'degrees_north'
    #                    long_name = 'latitude'
    lat = np.array(ncfile.variables['latitude'])

    # longitude
    #        Size:       1440x1
    #        Dimensions: longitude
    #        Datatype:   single
    #        Attributes:
    #                    units     = 'degrees_east'
    #                    long_name = 'longitude'
    lon = np.array(ncfile.variables['longitude'])

    # time
    #        Size:       nTimes x 1
    #        Dimensions: time
    #        Datatype:   int32
    #        Attributes:
    #                    units     = 'hours since 1900-01-01 00:00:00.0'
    #                    long_name = 'time'
    #                    calendar  = 'gregorian'
    timeDays = np.array(ncfile.variables['time'])/24

    ncfile.close()

    # Let's calculate the dateInt string from the filename
    # (note that these filenames were
    # set by hand after downloading from the Copernicus Climate Data Store)
    dateInt = np.empty_like(timeDays)
    YrDayStr = (fnFull.split('.')[-2]).split('-')[0]
    Yr0 = int(YrDayStr[0:4])
    Month0 = int(YrDayStr[4:])
    d0 = datetime.date(year=Yr0, month=Month0, day=1)
    extra_days = 0
    for ii in range(len(timeDays)):
        dateIntCur = d0 + datetime.timedelta(
            days=int(timeDays[ii]-timeDays[0])+extra_days)
        if calendar.isleap(dateIntCur.year) and dateIntCur.month == 2:
            extra_days = extra_days + 1

        # in YYYYMMDD format
        dateInt[ii] = int(dateIntCur.strftime('%Y%m') + '01')

    # net flux is sum of ssr and str (due to sign conventions),
    # and divide by seconds/day to convert from J/m^2 to W/m^2
    # because ssr and str are accumulations scaled to be over one day
    flux = (ssr + str) / 86400

    return flux, dateInt, lat, lon, timeDays


def read_ORAS5_sst(fnFull: str,
                   VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                           npt.NDArray, npt.NDArray,
                                           npt.NDArray]:
    """Read sea surface temperature values.

    Input: ORAS5 oceanic forcing reanalysis data in NetCDF format.

    https://www.ecmwf.int/en/forecasts/dataset/ocean-reanalysis-system-5
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # sosstsst
    #        Size:       360x180x1
    #        Dimensions: lon,lat,time_counter
    #        Datatype:   single
    #        Attributes:
    #                    standard_name      = 'Sea Surface temperature'
    #                    long_name          = 'Sea Surface temperature'
    #                    units              = 'C'
    #                    _FillValue         = 9.969209968386869e+36
    #                    missing_value      = 9.969209968386869e+36
    #                    online_operation   = 'ave(x)'
    #                    interval_operation = 1200
    #                    interval_write     = 2678400
    #                    offline_operation  = 'ave(x)'
    sst = np.array(ncfile.variables['sosstsst'])

    # lat
    #        Size:       180x1
    #        Dimensions: lat
    #        Datatype:   double
    #        Attributes:
    #                    standard_name = 'latitude'
    #                    long_name     = 'latitude'
    #                    units         = 'degrees_north'
    #                    axis          = 'Y'
    lat = np.array(ncfile.variables['lat'])

    # lon
    #        Size:       360x1
    #        Dimensions: lon
    #        Datatype:   double
    #        Attributes:
    #                    standard_name = 'longitude'
    #                    long_name     = 'longitude'
    #                    units         = 'degrees_east'
    #                    axis          = 'X'
    lon = np.array(ncfile.variables['lon'])

    # time_counter
    #        Size:       1x1
    #        Dimensions: time_counter
    #        Datatype:   double
    #        Attributes:
    #        standard_name = 'time'
    #        units         = 'seconds since 1979-01-16 00:00:00 UTC'
    #        calendar      = 'gregorian'
    # Convert from second to days
    timeDays = np.array(ncfile.variables['time_counter'])/86400

    ncfile.close()

    # Let's calculate the dateInt string from the filename
    dateInt = [int(fnFull.split('_')[-2] + '01')]

    return sst, dateInt, lat, lon, timeDays


def read_PSC_thickness(fnFull: str,
                       VERBOSE: bool) -> Tuple[npt.NDArray, npt.NDArray,
                                               npt.NDArray, npt.NDArray,
                                               npt.NDArray]:
    """Read sea ice thickness values from the Polar Science Center.
    http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/

    The binary thickness files were converted to NetCDF
    using Robbie Mallett's piomas_bin_reader:
    https://github.com/robbiemallett/piomas_bin_reader
    """
    # open NetCDF file
    ncfile = Dataset(fnFull)

    # provides some info on the variables, dimensions, etc.
    if VERBOSE:
        print(ncfile)

    # thickness
    #        Size:       120x360x12
    #        Dimensions: y,x,t
    #        Datatype:   double
    #        Attributes:
    #                    _FillValue  = NaN
    #                    coordinates = 'month lat lon'
    thickness = np.array(ncfile.variables['thickness'])

    # lat
    #        Size:       120x360
    #        Dimensions: y,x
    #        Datatype:   double
    #        Attributes:
    #                    _FillValue = NaN
    lat = np.array(ncfile.variables['lat'])

    # lon
    #        Size:       120x360
    #        Dimensions: y,x
    #        Datatype:   double
    #        Attributes:
    #                    _FillValue = NaN
    lon = np.array(ncfile.variables['lon'])

    # month
    #        Size:       12x1
    #        Dimensions: t
    #        Datatype:   int32
    month = np.array(ncfile.variables['month'])

    ncfile.close()

    # Let's calculate the dateInt string from the filename
    curYear = os.path.splitext(os.path.basename(fnFull))[0]
    dateInt = []
    timeDays = []
    for curMonth in month:
        dateInt.append(int(curYear + str(curMonth).zfill(2) + '01'))
        timeDays.append(0)

    return thickness, dateInt, lat, lon, timeDays
