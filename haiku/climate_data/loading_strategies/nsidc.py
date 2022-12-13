
from haiku.climate_data.loading_strategies.netcdf4_loading_strategy import \
    NetCDF4LoadingStrategy

from datetime import datetime
import numpy as np
import numpy.typing as npt
from netCDF4 import Dataset


class Era5T2mStrategy(NetCDF4LoadingStrategy):
    VARIABLE_KEY = "t2m"
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 1
    LONGITUDE_KEY_INDEX = 2

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)

        # cutoff temperature placeholders
        data[data > 1e20] = -100.0

        return data


class Oras5SstStrategy(NetCDF4LoadingStrategy):
    VARIABLE_KEY = "sosstsst"
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 1
    LONGITUDE_KEY_INDEX = 2

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)

        # cutoff temperature placeholders for plotting clarity
        data[data > 1e20] = -100.0

        return data


class PscThicknessStrategy(NetCDF4LoadingStrategy):
    VARIABLE_KEY = "thickness"
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 2
    LONGITUDE_KEY_INDEX = 1

    def _load_time(self, dataset: Dataset, key: str) -> npt.NDArray:
        # each timestep within this dataset file is a month
        months = self._load_variable_as_array(dataset, key)

        # this dataset has a string description which includes the year
        # so parse it here
        year = int(dataset.year.split(" ")[-1])

        # since this dataset does not have days
        # assume first day of the month for consistency
        day = 1
        dates = [datetime(year, month, day) for month in months]
        return np.asarray(dates)


class SicStrategy(NetCDF4LoadingStrategy):
    """Abstract. Not a valid loading strategy.

    Specific SIC strategies should implement this one.
    """
    VARIABLE_KEY = ""
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 1
    LONGITUDE_KEY_INDEX = 2

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        # Dimensions: (time, lat, lon)
        # Datatype: int16
        data = super()._load_data(dataset, key)

        # convert data to percent
        # and set any NaN values to 0.0
        data = 100.0 * data
        data[np.isnan(data)] = 0.0

        # round to nearest percent
        #TODO: enable rounding with: data = np.around(data)
        #These upper and lower thresholds don't currently do anything
        data[data < 0.0] = 0.0
        #sets land area to -100 to make auto plots easier to visualize
        data[data > 2500.0] = -100.0
        #These upper and lower thresholds don't currently do anything
        data[data > 100.0] = 100.0

        return data

    def _load_latitude(self, dataset: Dataset, key: str) -> npt.NDArray:
        # TODO: We should not have to hardcode here
        #       Requires a slight redesign to fix
        #       Should revisit if there is time
        try:
            return super()._load_latitude(dataset, key)
        except Exception:
            # this is a gross way to do this
            # but is significantly less code than creating polar
            # strategies for each one of these
            return super()._load_latitude(dataset, "ygrid")

    def _load_longitude(self, dataset: Dataset, key: str) -> npt.NDArray:
        # TODO: We should not have to hardcode here
        #       Requires a slight redesign to fix
        #       Should revisit if there is time
        try:
            return super()._load_longitude(dataset, key)
        except Exception:
            return super()._load_longitude(dataset, "xgrid")


class SicMergedStrategy(SicStrategy):
    VARIABLE_KEY = "cdr_seaice_conc_monthly"

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)


        return data


class SicNtStrategy(SicStrategy):
    VARIABLE_KEY = "nsidc_nt_seaice_conc_monthly"


class SicBtStrategy(SicStrategy):
    VARIABLE_KEY = "nsidc_bt_seaice_conc_monthly"
