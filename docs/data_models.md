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

Currently from the Monthly CESM Large Ensemble Community Project<sup>[KA15]</sup>, we extract several variables that we expect to be causally linked to Sea Ice concentration in the Arctic: 

   * **ICEFRAC**: the Sea Ice coverage fraction over each spatial grid element
   * **TS**: the atmospheric surface temperature 
   * **SST**: the sea surface temperature 
   * **HI**: the sea ice thickness
   
We have also begun including forcing terms into the model. 
These measurements and their impact on other observables are learned when the Koopman models are trained, but then when the KoopmanModel is used to forecast forward, a predefined time-series of the forcing terms can be passed in.  
Currently CO2 is treated as a forcing term and we use historically measured forcing and predicted CO2 based on anthropogenic mechanisms. This allows analysis of how reducing C02 can affect the FKPM predictions for Sea Ice concentration in the Arctic. 
A preliminary analysis can be found in the [initial results](../Results/initial_motivation_results) section.  More detail of the implementation of forcing terms in the Koopman model, can be found in the [Koopman Modeling](../koopman) section.


&nbsp;  
##Observational Data - NSIDC

Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS Passive Microwave Data, Version 1 | National Snow and Ice Data Center (nsidc.org)

We compare the results of the simulations with observational sea ice concentration data obtained from the National Snow and Ice Data Center. These data include gridded daily (every other day for SMMR data) and monthly averaged sea ice concentrations for both the north and south polar regions. The dataset includes data since 26 October 1978. The data shows the brightness temperature data derived from a few different sensors (microwave radiometers that sense emitted microwave radiation) that represent the sea ice concentration. The data are provided in the polar stereographic projection at a grid cell size of 25 x 25 km. For grid cell sizes of 25 x 25km, the pixel intensities represent the fractional amount of sea ice covering that cell (scaled by 250).

<figure>
<img src="../figs/results/NSIDC_data.png" alt="NSIDC data" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 2:</b> Sea ice concentrations in the northern and southern hemisphere from satellite data obtained from NSIDC. The data is shown in a polar projection and 250 on the color scale corresponds to 100% coverage.</figcaption>
</figure>  

Currently from the Monthly NSIDC datasets, we extract several variables that we expect to be causally linked to Sea Ice concentration in the Arctic: 

   * **seaice_conc**: the Sea Ice coverage fraction over each spatial grid element (NSIDC G02202_V3)<sup>[WA19]</sup>
   * **T2M**: Average Air Temperature at each spatial grid point 2 meters above sea surface (ERA5)<sup>[He19]</sup>
   * **SST**: Sea Surface Temperature: ocean temperature at surface of water (ORAS5)<sup>[Zu19]</sup>
   * **Sea Ice Thickness**: Sea Ice Thickness (PIOMAS v2.1)<sup>[SC11, ZH03]</sup>

&nbsp;  
##Cloud Cover Data - to support Value of New Data Estimator

In Phase II of the ACTM program, we aim to determine the value of potential new measurements that will have an outsized impact on improving the model forecasting accuracy, identification and characterization of tipping points, or improve model robustness.
If we are able to quantify the improvement to the HAIKU models of sea ice, particularly related to tipping points, we can then weigh the cost of collecting more measurements or developing better methods of estimating those measurements from current data.
Specifically, cloud cover measurements, under ice water temperature measurements, and atmospheric heat flux measurements are poorly estimated and rarely directly measured for the arctic region while likely having a large impact on the dynamics of sea ice concentration.
We will focus on cloud cover measurements for their potential to improve the sea ice concentration models in Phase II of DARPA ACTM.  
Further motivation of this choice and the associated analyses can be found in the [Value of New Data Estimator](https://bae-systems-haiku.github.io/HAIKU/analyses/#value-of-new-data-estimator-vonde) section of this document.

Considering the pros and cons of each dataset, we will primarily use cloud and surface heat flux data from ERA5 (1979-2022), MERRA-2 (1980-2022), Cloud-Aerosol Lidar and Infrared Pathfinder Satellite Observation (CALIPSO, 2006-2016), and Clouds and Earth's Radiant Energy System Energy Balanced And Filled (CERES-EBAF Ed4.0, 2006-2016) in ACTM Phase II to include cloud radiative influences in shaping the atmosphere-sea ice connection. ERA5 and MERRA-2 are considered the two best reanalysis products incorporating all available satellite and in-situ information and using the most updated 4- and 3-dimension assimilation schemes, respectively. 

Collecting accurate measurements of cloud cover are challenging and costly:

   * Most cloud cover measurements are derived from satelite imagery which is insufficient to estimate the elevation of the clouds accurately. 
   * Ship based measurements are far more accurate, but are incredibly costly to procure. 
   * Station data can be used as a cross-comparison, but is limited to regions on land or at least very near a station.
   * Additionally, procuring measurements can be difficult within the zone of influence of certain countries. 
   
By accurately assessing the value of new measurements with specific locations, temporal frequencies, and specified precision can provide the needed context to properly balance the direct cost of making the measurements with the improvement to climate forecasting of Sea Ice tipping points in the Arctic region.

##Dataset Access

Instructions to download and preprocess datasets currently in use on HAIKU can be found on the [HAIKU github](https://bae-systems-haiku.github.io/HAIKU/data_models/). 
This page will continue to be updated as the HAIKU code is built out, validated, and approved for public use and will contain instructions to get HAIKU operational on your system.


##Data Citations

#####CESM Large Ensemble data

<sup>[KA15]</sup>  Kay, J. E., Deser, C., Phillips, A., Mai, A., Hannay, C., Strand, G., Arblaster, J., Bates, S., Danabasoglu, G., Edwards, J., Holland, M. Kushner, P., Lamarque, J.-F., Lawrence, D., Lindsay, K., Middleton, A., Munoz, E., Neale, R., Oleson, K., Polvani, L., and M. Vertenstein (2015), The Community Earth System Model (CESM) Large Ensemble Project: A Community Resource for Studying Climate Change in the Presence of Internal Climate Variability, Bulletin of the American Meteorological Society, doi: 10.1175/BAMS-D-13-00255.1, 96, 1333-1349

#####NSIDC Sea Ice Concentration data

<sup>[WA19]</sup> Walsh, J. E., W. L. Chapman, F. Fetterer, and J. S. Stewart. 2019. Gridded Monthly Sea Ice Extent and Concentration, 1850 Onward, Version 2. [NSIDC G02202_V3]. Boulder, Colorado USA. NSIDC: National Snow and Ice Data Center. doi: https://doi.org/10.7265/jj4s-tq79. [Date Accessed].

#####ERA5 Monthly

<sup>[HE19]</sup> Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I., Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2019): ERA5 monthly averaged data on single levels from 1979 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). (Accessed on < DD-MMM-YYYY >), 10.24381/cds.f17050d7

#####ORAS5

<sup>[ZU19]</sup> Zuo, H., Balmaseda, M. A., Tietsche, S., Mogensen, K., and Mayer, M.: The ECMWF operational ensemble reanalysis–analysis system for ocean and sea ice: a description of the system and assessment, Ocean Sci., 15, 779–808, https://doi.org/10.5194/os-15-779-2019, 2019.

#####PIOMAS Volume time series and uncertainties

<sup>[SC11]</sup> Schweiger, A., R. Lindsay, J. Zhang, M. Steele, H. Stern, Uncertainty in modeled arctic sea ice volume, J. Geophys. Res., doi:10.1029/2011JC007084, 2011

#####PIOMAS Model details

<sup>[ZH03]</sup> Zhang, J.L. and D.A. Rothrock, “Modeling global sea ice with a thickness and enthalpy distribution model in generalized curvilinear coordinates“, Mon. Weather Rev., 131, 845-861, 2003

