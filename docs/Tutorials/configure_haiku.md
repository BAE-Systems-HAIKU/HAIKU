
#Configuring HAIKU on your System

##Cloning HAIKU
First, clone the github HAIKU repository:

```
git clone https://github.com/BAE-Systems-HAIKU/HAIKU.git HAIKU
cd HAIKU
```

##Setting up python env
Next setup the python environment:
HAIKU expects a Linux environment running Python 3.8 or higher. We recommend using a Python Virtual Environment to isolate HAIKU dependencies from the rest of your system. You can use your favorite approach to set this up, but we provide instructions assuming python venv and pip as environment and package managers. There are Python library dependencies users need to download before running the system. To do so, users should use `pip` in conjunction with the `requirements.txt` file found in the base directory of the public GitHub repository. 
To create this environment and install dependencies, run:

```
#create the environment
python3.8 -m venv ./haiku-venv
#activate the python environment
source ./haiku-venv/bin/activate
#install HAIKU dependencies in local virtual environment
pip install -r requirements.txt
```

Finally, the system `PYTHONPATH` environment variable must include the root directory of the haiku software. For example, if this codebase was located at `/home/test/core/haiku`, it could be added to the system `PYTHONPATH` with the following command:
```
export PYTHONPATH=$PYTHONPATH:/home/test/core/haiku
```
This command can be added to `~/.bashrc` on a Linux system so that it is applied automatically whenever a terminal is launched.


## Installing CDO

We execute this data alignment with the Climate Data Operators (cd0) tool available [here](https://code.mpimet.mpg.de/projects/cdo/). This is required to compare observational data to CESM data, this is a default part of the HAIKU analyses.  This code must be compiled before it can be run and requires the c++14 standard. Additionally, the data is in netCDF format, which means cdo must be compiled with [netDCF](https://github.com/Unidata/netcdf-c) which itself requires [hdf5](https://github.com/HDFGroup/hdf5).  To accomplish the remapping from stereographic projection cdo requires [proj](https://proj.org/) (which requires a new enough version of [sqlite](https://github.com/sqlite) to compile.  `make_install_cdo.sh` should handle downloading and compiling all the appropriate files as well as executing the appropriate regridding of all netCDF files. CDO requires an input parameter file to specific details of the target grid.  The required target parameter file is given in `configs/target_grid.txt`. The NSIDC sea ice concentration data V4 has an underspecified grid defined in its .nc data file, so we also include the fully defined northern hemisphere and southern hemisphere grid definition for the stereographic grid (`configs/nh_icefrac_grid.txt`) and southern hemisphere (`configs/sh_icefrac_grid.txt`). The `data_downloader.sh` script will automatically regrid the NSIDC data.

On Ubuntu 20.04, this binary can be downloaded directly from the package repositories with `apt install cdo`. If you prefer, the binary can be downloaded directly from the website linked above. If choosing this route, ensure the downloaded binary is placed on your system path (e.g., in `/usr/bin`). HAIKU expects a globally accesible CDO binary.


## Data Download

The data currently used in this system is described in detail in [HAIKU data](../../data_models) . We include a `data_downloader.sh` that automates downloading of all the CESM1 (simulation) data and some of the observation/reanalysis data. Unfortunately, much of the reanalysis data must be downloaded by hand, but descriptions of exactly what parameters to use to download the data are listed in the `data_downloader.sh` script. We will release a set of preprocessed data in the near future to reduce user effort and avoid possible issues in configuring the data, but will also maintain description of the full process for reproduceability. 

### NCAR data
We extract several variables from the NCAR getway hosted [here](https://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html).
This data is free to download, but required registering a user account first. With a user account, an API Token
is provided to access data hosted by NCAR/earthsystemsgrid.org. You can access this token by browsing to
"Account Home" (listed under your username). This API token must be provided to `data_downloader.sh` as the single argument.

### NSIDC data
The NSIDC data are hosted [here](ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_V4/) covering 1987-present
at a monthly temporal resolution. This data is hosted on an ftp server and no additional login information
is required to download it. However, the coordinate system of the NSIDC data do not align with the CESM1 data. Currently, the preprocessing step will convert coordinate systems using python CDO bindings to (default) polar stereographic or lon-lat coordinates. These steps are described in [Running HAIKU](../running_the_code#data-preprocessing-and-climatedata-objects).
