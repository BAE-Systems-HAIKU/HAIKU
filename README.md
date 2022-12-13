# Initial setup
## Git initialization

### Preferred initial setup with ssh public key in bitbucket account:

`git clone ssh://git@git.dmz.devlnk.net/darpaactm/core.git HAIKU`

### Alternative approach using https:

`git clone https://git.dmz.devlnk.net/scm/darpaactm/core.git HAIKU`

*note, if this fails, you may need to temporarily disable the git
 sslVerification as below:*

`git config --global http.sslVerify "false"`

*and also ignore a corporate proxy when connecting to the sdmz's bitbucket
 repo:*

`export no_proxy=git.dmz.devlnk.net`

*(the above line can likely go in your ~/.bashrc file so you don't need to run it
 every time.)*

``$cd HAIKU``

# Set up Python env

Set up a python virtual environments with the supplied `requirements.txt` file (or use whatever mechanism you prefer to set up your working python environment).

```
python3 -m venv ./venv
source ./venv/bin/activate
pip install -r requirements.txt
```

This virtual environment can now be activated at any time, and the installed libraries are kept separate
from the global libraries.

# Other Dependencies

`cdo` is required for running this software. On Ubuntu, the binaries can be easily installed with: `apt install cdo`

Specific releases can be downloaded from the website: https://code.mpimet.mpg.de/projects/cdo/files

Building CDO on CentOS 7:

```
yum install epel-release
yum install netcdf netcdf-devel
wget https://code.mpimet.mpg.de/attachments/download/23323/cdo-1.9.9.tar.gz
tar xzf cdo-1.9.9.tar.gz
cd cdo-1.9.9
./configure --with-netcdf
make
make install
```

# HAIKU software structure

 - **climate**: contains code to run or interface directly with climate
   data (including preprocessing and filtering).

 - **koopman**: contains code to generate/train Koopman operator based models

 - **plotting**: contains code to generate plots derived from koopman or
   climate data objects

 - **prediction**: contains code to make predictions with koopman models

 - **analytics**: contains code to process models or models outputs to
   generate various analyses.

 - **utilities**: methods used across the other systems, typically related to
   data transformations and string parsing; also includes configuration class
   and defaults

Sea Ice Data Downloading
========================
The data currently used in this system is described in detail in the full documentation
linked below. We include a `data_downloader.sh` that automates downloading of all the CESM1 (simulation)
data and some of the observation/reanalysis data. Unfortunately, much of the reanalysis data must be
downloaded by hand, but descriptions of exactly what parameters to use to download the data are listed in the
`data_downloader.sh` script.

Requirements
------------
In order for `data_downloader.sh` to run, you must have wget installed and a compiler that is c++14 compatible.
This script has been tested under RedHat Linux 7.9 with gcc 8.3 and Ubuntu 20.04. Versions and all dependencies
were selected such that they can all be built with make and other build systems should not be required.

NCAR data
---------
We extract several variables from the NCAR getway hosted [here](https://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html).
This data is free to download, but required registering a user account first. With a user account, an API Token
is provided to access data hosted by NCAR/earthsystemsgrid.org. You can access this token by browsing to
"Account Home" (listed under your username). This API token must be provided to `data_downloader.sh` as the single argument.

NSIDC data
----------
The NSIDC data are hosted [here](ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_V4/) covering 1987-present
at a monthly temporal resolution. This data is hosted on an ftp server and no additional login information
is required to download it. However, the coordinate system of the NSIDC data do not align with the CESM1 data
and currently must be converted to latitude/longitude to make a direct comparison to the CESM1 output.

#### Regridding data with Climate Data Operators
We execute this data alignment with the Climate Data Operators (cd0) tool available [here](https://code.mpimet.mpg.de/projects/cdo/).
This code must be compiled before it can be run and requires the c++14 standard. Additionally, the data
is in netCDF format, which means cdo must be compiled with [netDCF](https://github.com/Unidata/netcdf-c) which itself requires [hdf5](https://github.com/HDFGroup/hdf5).
To accomplish the remapping from stereographic projection cdo requires [proj](https://proj.org/) (which
requires a new enough version of [sqlite](https://github.com/sqlite) to compile.
`data_downloader.sh` should handle downloading and compiling all the appropriate files as well as executing
the appropriate regridding of all netCDF files. CDO requires an input parameter file to specific details of the target grid.
The required target parameter file is given in `configs/target_grid.txt`. The NSIDC sea ice concentration data V4 has an underspecified
grid defined in its .nc data file, so we also include the fully defined northern hemisphere and southern hemisphere grid definition
for the stereographic grid (`configs/nh_icefrac_grid.txt`) and southern hemisphere (`configs/sh_icefrac_grid.txt`).

#### Quickstart data download and preprocess:

   1. Navigate to [NCAR User API token](https://www.earthsystemgrid.org/ac/user/apiTokenDisplay.html) to copy your API token (you will need to make and account/login)

   2. `bash data_downloader.sh pasted_api_token` to download and preprocess the relevant data.

   3. The resultant files are split into the CESM1/ and regridded_nsidc/ directories which are ready to be ingested by the HAIKU code

   4. `data_downloader.sh` will also spit out text describing urls and rough instructions for any additional observational/reanalysis data that we've used in our analysis

# Full Documentation for HAIKU

Please find the full documentation of this repository along with recent results and future work
at [github.io](bae-systems-haiku.github.io/HAIKU/)