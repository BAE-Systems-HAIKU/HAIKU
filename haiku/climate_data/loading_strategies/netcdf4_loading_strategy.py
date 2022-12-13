from __future__ import annotations

import logging

import numpy as np
import numpy.typing as npt
from cftime import num2date
from netCDF4 import Dataset, Variable

from haiku.climate_data.climate_data import ClimateData


class NetCDF4LoadingStrategy:
    """Base class for loading strategies to implement.

    A NetCDF4 loading strategy extracts a single variable from a
    NetCDF4 (*.nc) dataset file.

    Along with the variable data, these strategies also extract
    latitude, longitude, and dates for each timestep.

    Since NetCDF4 files can contain multiple variables, sometimes
    multiple strategies will need to be used in order to extract
    all relevant data for a given use case.

    Just call load with a path to a dataset and this class should
    handle the rest.
    """
    VARIABLE_KEY = ""
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 0
    LONGITUDE_KEY_INDEX = 0

    def __eq__(self, o: NetCDF4LoadingStrategy) -> bool:
        return self.VARIABLE_KEY == o.VARIABLE_KEY

    def __hash__(self):
        return hash(self.VARIABLE_KEY)

    def _load_variable(self, dataset: Dataset,
                       key: str) -> Variable:
        try:
            return dataset.variables[key]
        except Exception as e:
            logging.info("Valid variables: %s", list(dataset.variables))
            raise e

    def _load_variable_as_array(self, dataset: Dataset,
                                key: str) -> npt.NDArray:
        return np.array(self._load_variable(dataset, key))

    def _load_latitude(self, dataset: Dataset, key: str) -> npt.NDArray:
        return self._load_variable_as_array(dataset, key)

    def _load_longitude(self, dataset: Dataset, key: str) -> npt.NDArray:
        return self._load_variable_as_array(dataset, key)

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        """Override this method if additional processing is needed."""
        # TODO: Do we want to reverse the latitudes here?
        return self._load_variable_as_array(dataset, key)

    def _load_time(self, dataset: Dataset, key: str) -> npt.NDArray:
        timesteps = self._load_variable(dataset, key)
        time_units = timesteps.units
        time_calendar = timesteps.calendar

        # here we convert the time variable to python datetime
        # if python datetime is not possible, it will be converted
        # to a cftime object instead
        time_values = np.array(timesteps)
        real_time_values = num2date(
            time_values,
            time_units,
            calendar=time_calendar,
            only_use_cftime_datetimes=False
        )
        return real_time_values

    def _extract_data(self, dataset: Dataset) -> ClimateData:
        # extract all info from dataset
        data = self._load_data(dataset, self.VARIABLE_KEY)

        data_dimensions = dataset.variables[self.VARIABLE_KEY].dimensions

        latitude_key = data_dimensions[self.LATITUDE_KEY_INDEX]
        latitudes = self._load_latitude(dataset, latitude_key)

        longitude_key = data_dimensions[self.LONGITUDE_KEY_INDEX]
        longitudes = self._load_longitude(dataset, longitude_key)

        timestep_key = data_dimensions[self.TIMESTEP_KEY_INDEX]
        timesteps = self._load_time(dataset, timestep_key)

        # build climate data object and store in list
        climate_data = ClimateData(
            data,
            timesteps,
            latitudes,
            longitudes,
            description=self.VARIABLE_KEY)

        return climate_data

    def load_from_disk(self, filepath: str) -> ClimateData:
        with Dataset(filepath) as dataset:
            return self._extract_data(dataset)

    def load_from_memory(self, dataset: Dataset) -> ClimateData:
        return self._extract_data(dataset)
