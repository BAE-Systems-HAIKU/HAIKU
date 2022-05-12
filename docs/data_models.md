# Data and Climate models
The HAIKU program will focus on analyzing decade-scale predictions of Arctic sea ice concentrations. 
We plan to initially use the pre-generated data across a variety of CICE4 parameters to initially construct the HAIKU models, but expect to need to run the CICE4 model eventually before transitioning to the full CESM model. 
We also plan to leverage real measurement data to train our Hybrid Koopman-Climate Model which will learn to model dynamics present in real data that is not present in the models and their associated pre-generated data. 

To enable rapid early results and validate our approach, we begin with a stand-alone sea ice model: the Los Alamos sea ice model [CICE4](https://www.cesm.ucar.edu/models/cesm1.0/cice/) and associated [documentation](http://www.ccpo.odu.edu/~klinck/Reprints/PDF/cicedoc2015.pdf).

At the midpoint of Phase 1, we plan to move to the full Community Earth System Model (CESM2) to better model the coupled climate variability. CESM2 integrates CICE5 along with atmosphere (CAM6), ocean (POP2), land (CLM5), and ice sheet (CISM2) models of the NCAR modeling framework.

##Initial Climate Models and Datasets

The CICE5 model requires atmospheric data, including: monthly downward shortwave radiation data at the surface, precipitation, cloud fraction, and four times daily data of 2m air temperature, 2m specific humidity, 10m wind, and 10m air density. 
Most of the atmospheric data will be derived from the ERA5 reanalysis accessible via the ECWMF data portal ([ERA5](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5)).

The CICE5 model also requires oceanic data, including: monthly sea surface temperature, sea surface salinity, sea surface height, mixed layer depth, and ocean currents. 
The oceanic forcing will come from several indicators available at the ECMWF Ocean Reanalysis System 5 ([ORAS5 via ECWMF](https://www.ecmwf.int/en/forecasts/dataset/ocean-reanalysis-system-5)).
In addition, the CICE5 model requires data for vertical heat transport from the deep ocean, obtainable as CESM model output available via [NCAR climate data gateway](https://www.earthsystemgrid.org/). 

The first step is to confirm sensitivity observations concluded by the climatology community. 
Before investigative sensitivity, we will first need to find a correspondence between model outputs and observational data for validation. Observational data for sea ice concentration can be obtained from the [National Snow and Ice Data Center](https://nsidc.org/data/NSIDC-0051/versions/1)  (NSIDC).
We will need to work in a common projection and resolution to avoid interpolation.


##Model Prediction - CESM Large Ensemble Project
The models are computationally expensive to run. Instead, we look at the CESM Large Ensemble Community Project. This project has produced a publicly available set of climate model simulations performed with the nominal 1-degree latitude/longitude version of the Community Earth System Model version 1 (CESM1) with CAM5.2 as its atmospheric component [CESM data](https://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html). 

<figure>
<img src="../figs/results/CESM_data.png" alt="CESM data" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 1:</b> Sea ice concentration percentage in the northern hemisphere generated from the CESM Large Ensemble Community Project. The data is shown in longitude-latitude coordinates. </figcaption>
</figure>

##Observational Data - NSIDC

Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS Passive Microwave Data, Version 1 | National Snow and Ice Data Center (nsidc.org)

We compare the results of the simulations with observational sea ice concentration data obtained from the National Snow and Ice Data Center. These data include gridded daily (every other day for SMMR data) and monthly averaged sea ice concentrations for both the north and south polar regions. The dataset includes data since 26 October 1978. The data shows the brightness temperature data derived from a few different sensors (microwave radiometers that sense emitted microwave radiation) that represent the sea ice concentration. The data are provided in the polar stereographic projection at a grid cell size of 25 x 25 km. For grid cell sizes of 25 x 25km, the pixel intensities represent the fractional amount of sea ice covering that cell (scaled by 250).

<figure>
<img src="../figs/results/NSIDC_data.png" alt="NSIDC data" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 2:</b> Sea ice concentrations in the northern and southern hemisphere from satellite data obtained from NSIDC. The data is shown in a polar projection and 250 on the color scale corresponds to 100% coverage.</figcaption>
</figure>


##Dataset Access

Instructions to download and preprocess datasets currently in use on HAIKU can be found on the [HAIKU github](https://github.com/BAE-Systems-HAIKU/HAIKU). 
This page will continue to be updated as the HAIKU code is built out, validated, and approved for public use and will contain instructions to get HAIKU operational on your system.




