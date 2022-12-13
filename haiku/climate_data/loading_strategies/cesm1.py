
import numpy as np
import numpy.typing as npt
from haiku.climate_data.loading_strategies.netcdf4_loading_strategy import \
    NetCDF4LoadingStrategy
from netCDF4 import Dataset


class IceFracStrategy(NetCDF4LoadingStrategy):
    """Sea Ice Concentration from CESM1 LENS."""
    VARIABLE_KEY = "ICEFRAC"
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 1
    LONGITUDE_KEY_INDEX = 2

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)

        # convert data to percent
        data = 100.0 * data

        #TODO: enable rounding with: data = np.around(data)

        # ensure percentage is between 0 and 100
        data[data > 100.0] = -100.0
        data[data < 0.0] = 0.0

        return data

    def _load_latitude(self, dataset: Dataset, key: str) -> npt.NDArray:
        return super()._load_latitude(dataset, key)


class SstStrategy(NetCDF4LoadingStrategy):
    """Sea Surface Temperature from CESM1 LENS."""
    VARIABLE_KEY = "SST"
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 2
    LONGITUDE_KEY_INDEX = 3

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)

        # permute dimensions to match icefrac dimensions
        data = np.squeeze(np.transpose(data, (0, 2, 3, 1)))

        # cutoff temperature placeholders
        data[data > 1e20] = -100.0

        return data


class TsStrategy(NetCDF4LoadingStrategy):
    """Atmospheric Surface Temperature from CESM1 LENS."""
    VARIABLE_KEY = "TS"
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 1
    LONGITUDE_KEY_INDEX = 2

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)

        # cutoff temperature placeholders
        data[data > 1e20] = 0.0

        return data


class HiStrategy(NetCDF4LoadingStrategy):
    """Ice Volume per Unit Grid Cell Area from CESM1 LENS."""
    VARIABLE_KEY = "hi"
    TIMESTEP_KEY_INDEX = 0
    LATITUDE_KEY_INDEX = 1
    LONGITUDE_KEY_INDEX = 2

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)

        # cutoff value placeholders
        data[data > 1e20] = -100.0

        return data


class ForcingStrategy(NetCDF4LoadingStrategy):
    """Abstract. Implement by specific forcing terms."""

    def _load_data(self, dataset: Dataset, key: str) -> npt.NDArray:
        data = super()._load_data(dataset, key)

        data = np.expand_dims(data, axis=(1, 2))

        return data

    def _load_latitude(self, dataset: Dataset, key: str) -> npt.NDArray:
        return super()._load_latitude(dataset, "lat")

    def _load_longitude(self, dataset: Dataset, key: str) -> npt.NDArray:
        return super()._load_longitude(dataset, "lon")


class Ch4vmrStrategy(ForcingStrategy):
    """Greenhouse Gas Forcing from CESM1 LENS."""
    VARIABLE_KEY = "ch4vmr"
    TIMESTEP_KEY_INDEX = 0


class Co2vmrStrategy(ForcingStrategy):
    """Greenhouse Gas Forcing from CESM1 LENS."""
    VARIABLE_KEY = "co2vmr"
    TIMESTEP_KEY_INDEX = 0


class F11vmrStrategy(ForcingStrategy):
    """Greenhouse Gas Forcing from CESM1 LENS."""
    VARIABLE_KEY = "f11vmr"
    TIMESTEP_KEY_INDEX = 0


class F12vmrStrategy(ForcingStrategy):
    """Greenhouse Gas Forcing from CESM1 LENS."""
    VARIABLE_KEY = "f12vmr"
    TIMESTEP_KEY_INDEX = 0


class N2ovmrStrategy(ForcingStrategy):
    """Greenhouse Gas Forcing from CESM1 LENS."""
    VARIABLE_KEY = "n2ovmr"
    TIMESTEP_KEY_INDEX = 0


class SolTsiStrategy(ForcingStrategy):
    """Greenhouse Gas Forcing from CESM1 LENS."""
    VARIABLE_KEY = "sol_tsi"
    TIMESTEP_KEY_INDEX = 0
