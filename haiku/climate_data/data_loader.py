

import logging
import os
from typing import List

from haiku.climate_data.climate_data import ClimateData
from haiku.climate_data.loading_strategies.cesm1 import (Ch4vmrStrategy,
                                                         Co2vmrStrategy,
                                                         F11vmrStrategy,
                                                         F12vmrStrategy,
                                                         HiStrategy,
                                                         IceFracStrategy,
                                                         N2ovmrStrategy,
                                                         SolTsiStrategy,
                                                         SstStrategy,
                                                         TsStrategy)
from haiku.climate_data.loading_strategies.netcdf4_loading_strategy import \
    NetCDF4LoadingStrategy
from haiku.climate_data.loading_strategies.nsidc import (Era5T2mStrategy,
                                                         Oras5SstStrategy,
                                                         PscThicknessStrategy,
                                                         SicBtStrategy,
                                                         SicMergedStrategy,
                                                         SicNtStrategy)
from netCDF4 import Dataset


class DataLoader:

    SUPPORTED_LOADING_STRATEGIES = [
        SstStrategy,
        TsStrategy,
        IceFracStrategy,
        HiStrategy,
        Ch4vmrStrategy,
        Co2vmrStrategy,
        F11vmrStrategy,
        F12vmrStrategy,
        N2ovmrStrategy,
        SolTsiStrategy,
        Era5T2mStrategy,
        Oras5SstStrategy,
        SicMergedStrategy,
        SicBtStrategy,
        SicNtStrategy,
        PscThicknessStrategy
    ]

    def __variable_exists(self, class_variable: str,
                          dataset_variables: List[str]) -> bool:
        return class_variable in dataset_variables

    def determine_strategies(self, filepath: str) \
            -> List[NetCDF4LoadingStrategy]:
        strategies = []
        # extract variable names from the dataset
        with Dataset(filepath) as dataset:
            variable_names = list(dataset.variables)

        # include any supported loading strategy which has its
        # variable key included in this dataset's variables
        for loading_strategy in self.SUPPORTED_LOADING_STRATEGIES:
            if self.__variable_exists(loading_strategy.VARIABLE_KEY,
                                      variable_names):
                strategies.append(loading_strategy())

        if len(strategies) == 0:
            err = f"Unsupported variables: {variable_names}"
            logging.error(err)
            raise ValueError(err)

        return strategies

    def load_directory(self, directory: str) -> List[ClimateData]:
        if not os.path.exists(directory):
            err = f"Directory not found: {directory}"
            logging.error(err)
            raise ValueError(err)

        filenames = os.listdir(directory)
        if len(filenames) == 0:
            err = f"No files found in directory: {directory}"
            logging.warning(err)
            return []

        results: List[ClimateData] = []
        for filename in os.listdir(directory):
            filepath = os.path.join(directory, filename)
            try:
                loaded_data = self.load_file(filepath)
                results.extend(loaded_data)
            except Exception:
                logging.warn("Unable to load. Skipping file: %s", filepath)

        return results

    def load_file(self, filepath: str) -> List[ClimateData]:
        if not os.path.exists(filepath):
            err = f"Filepath not found: {filepath}"
            logging.error(err)
            raise ValueError(err)

        load_strategies = self.determine_strategies(filepath)
        results = []
        for load_strategy in load_strategies:
            loaded_data = load_strategy.load(filepath)
            results.append(loaded_data)
        return results
