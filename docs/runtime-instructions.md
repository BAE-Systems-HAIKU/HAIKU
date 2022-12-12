# Runtime Instructions

## Overview

This document provides an overview of the HAIKU runtimes and supported data formats.

## Data Sources and Variables

The following are a list of currently supported data sources and variables.

__Observational (NSIDC)__
- Sea Ice Concentration
    - Goddard Merged
    - Goddard NASA Team
    - Goddard Bootstrap
- Sea Surface Temperature (ORAS5)
- 2 Meter Temperature (ERA5)
- Thickness (PIOMAS)

__Simulated (CESM1)__
- Sea Ice Concentration
- Sea Surface Temperature
- Surface Temperature
- Forcing
    - ch4 volume mixing ratio
    - co2 volume mixing ratio
    - f11 volume mixing ratio
    - f12 volume mixing ratio
    - n20 volume mixing ratio
    - Total solar irradiance (W/m2)

## Preprocessing

The first step of running HAIKU is data preprocessing.

In this step, HAIKU extracts variables from the data source datasets and formats it for training. The software accomplishes this with the following steps:

1. Loading
    - Extract variables of interest from dataset, load them into memory, and discard the rest

2. Sorting
    - Data is sorted chronologically once loaded into memory

3. Temporal Interpolation
    - If data is missing timesteps, the data is interpolated linearly and a log message is recorded
    - If data is on a different lat/lon grid

4. Aggregation
    - There are two supported runtimes for aggregate:
        - Combine
            - Data for the same variable types are combined into the same data structure
            - Use case: Some datasets are split into multiple files when stored on disk
            - Assumption: Data does not overlap timesteps across the split files
        - Average
            - All loaded data is averaged into a single dataset
            - Use case: CESM1 Large Ensemble datasets
            - Assumption: All data being averaged is the same variable and has the same number of timesteps

5. Spatial Interpolation
    - Data is remapped to match a prespecified grid

6. Serialization
    - Formatted data is serialized to disk as a `pickle` file
    - This output is used subsequently as input to the next step
