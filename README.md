Sea Ice Data Downloading
========================
The data used so far in this repository comes from two different sources NCAR and NSIDC. NCAR provides simulated data from multiple runs from the CESM1 model, whereas NSIDC provides observational data of sea ice concentration. The script `data_downloader.sh` aims to automate downloading and assembling this data as much as possible for use in fitting Koopman models. There are a few steps required to obtain this data.


Requirements
------------
In order for `data_downloader.sh` to run you must have wget installed and compiler that is c++14 compatible. This script has been tested under RedHat Linux 7.9 with gcc 8.3 and Ubuntu 20.04. Versions all dependencies were selected such that they can all be built with make and other build systems should not be required. 


NCAR data
---------
We extract only one variable from the full CESM set of variables this point in time: ICEFRAC this data is hosted [here](https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm4.cesmLE.atm.proc.monthly_ave.ICEFRAC.html). This data is free to download but requires registering a user account first. With a user account, an API Token is provided to access data hosted by NCAR/earthsystemgrid.org. You can access this token by browsing to "Account Home" (listed under your username). This API token must be provided to `data_downloader.sh` as the single argument:

NSIDC data
----------

The NSIDC data are hosted [here](ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_V4/) covering 1987-present at a monthly temporal resolution. This data is hosted on an ftp server and no additional login information is required to download this data. However, the coordinate systems of the NSIDC data do not align with the CESM1 data and must be converted to latitude/longitude to make a direct comparison to the CESM1 output.

#### Regridding data with Climate Data Operators
We execute this data alignment with the Climate Data Operators (cdo) tool available [here](https://code.mpimet.mpg.de/projects/cdo/). This code must be compiled before it can be run requires the c++14 standard. Additionally the data is in netCDF format, which means cdo must be compiled with [netCDF](https://github.com/Unidata/netcdf-c) which itself requires [hdf5](https://github.com/HDFGroup/hdf5). To accomplish the remapping from stereographic projection cdo requires [proj](https://proj.org/) (which requires a new enough version of [sqlite](https://github.com/sqlite) to compile. `data_downloader.sh` should handle downloading and compiling all the appropriate files as well as executing the appropriate regridding of all netCDF files. CDO requires an input parameter file to specific details of the target grid. The required parameter file is given in cesm1_grid.txt.

#### Miscellaneous
The NSIDC data has a gap of two months contained squarely in the middle of the period of interest. To handle this we linearly interpolated between the gaps and these files are uploaded separately as:

* latlon_seaice_conc_monthly_nh_zzz_198712_v03r01Interp.nc
* latlon_seaice_conc_monthly_nh_zzz_198712_v03r01Interp.nc
* latlon_seaice_conc_monthly_sh_zzz_198712_v03r01Interp.nc
* latlon_seaice_conc_monthly_sh_zzz_198801_v03r01Interp.nc

#### Quickstart

   1. Navigate to [NCAR User API token](https://www.earthsystemgrid.org/ac/user/apiTokenDisplay.html) to copy your API token (you will need to make an account/login).

   2. `bash data_downloader.sh paste_your_api_token` to download and preprocess the relevant data.

   3. The resultant files are split into the CESM1/ and regridded_nsidc/ directories which are ready to be ingested by the HAIKU code

As we continue processing and training HAIKU these instructions and scripts will be updated to include the additional variables of interest. 
