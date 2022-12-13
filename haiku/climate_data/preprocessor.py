"""
@classification
UNCLASSIFIED

@itar
ITAR CONTROLLED

@copyright
Copyright BAE Systems
Copyright (c) 2022 AIMdyn Inc.

Please reference the LICENSE.txt file included with this software
package for all terms and conditions related to this software
including license, rights and distribution.
"""

import logging
import os
from typing import List

import matplotlib.pyplot as plt
from haiku.cdo_wrapper.cdo_wrapper import CdoWrapper
from haiku.cdo_wrapper.grid_information import GridParser
from haiku.climate_data.climate_data import ClimateData
from haiku.climate_data.data_loader import DataLoader
from haiku.climate_data.interpolation import Interpolator
from haiku.climate_data.loading_strategies.netcdf4_loading_strategy import \
    NetCDF4LoadingStrategy
from haiku.climate_data.loading_strategies.cesm1 import ForcingStrategy
from haiku.climate_data.loading_strategies.nsidc import SicStrategy


class Preprocessor:

    def __init__(self, reference_grid_filepath: str):
        self.data_loader = DataLoader()
        self.interpolator = Interpolator()
        self.cdo_wrapper = CdoWrapper()
        self.grid_filepath = reference_grid_filepath
        self.grid_parser = GridParser()
        self.grid = self.grid_parser.parse_cdo_gridfile(self.grid_filepath)

    def _get_all_files_in_directory(self, directory: str) -> List[str]:
        # get all files from input directory
        filepaths = []
        for dirpath, _, filenames in os.walk(directory):
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                filepaths.append(filepath)
            # break because we want to ignore any subdirectories here
            break
        return filepaths

    def _determine_loading_strategies(self, filepaths: List[str]) \
            -> List[NetCDF4LoadingStrategy]:
        strategies = set()
        for filepath in filepaths.copy():
            try:
                current_strategies =\
                    self.data_loader.determine_strategies(filepath)
                for current_strategy in current_strategies:
                    strategies.add(current_strategy)
            except Exception as e:
                logging.warning("Unable to determine strategies from file. "
                                "Skipping: %s", filepath)
                filepaths.remove(filepath)
                logging.exception(e)
        return list(strategies)

    def _load(self, input_directory: str) -> List[ClimateData]:
        logging.info("Loading data from %s", input_directory)
        data = self.data_loader.load_directory(input_directory)

        if len(data) == 0:
            logging.warning("No data was loaded")
            return data

        logging.info("Loaded %d datafiles into memory", len(data))
        return data

    def _sort(self, data: List[ClimateData]) -> List[ClimateData]:
        logging.info("Sorting %d datapoints...", len(data))
        # sort all the data by the first date
        data.sort(key=lambda e: e.dates[0])
        logging.info("Sorted data spans from %s to %s",
                     data[0].dates[0], data[-1].dates[-1])
        return data

    def _temporal_interpolation(self,
                                data: List[ClimateData]) -> List[ClimateData]:
        logging.info("Starting temporal interpolation step...")
        # interpolate if necessary
        if len(data) > 1:
            interpolated_data =\
                self.interpolator.interpolate_climate_data(data)
            logging.info("Interpolated %d new datapoints",
                         len(interpolated_data))
        else:
            logging.info("Only 1 datapoint so no temporal interpolation")
            interpolated_data = []

        # if we have interpolated data, insert it and sort the list again
        if len(interpolated_data) > 0:
            data.extend(interpolated_data)
            data = self._sort(data)

        return data

    def _combine(self, data: List[ClimateData]) -> ClimateData:
        logging.info("Combining %d datapoints into one object", len(data))
        # combine the sorted ClimateData objects into one
        # containing all timesteps
        master_climate_data = ClimateData()
        for datapoint in data:
            master_climate_data.extend(datapoint)

        logging.info("Master Climate Data has %d timesteps",
                     master_climate_data.data.shape[0])

        return master_climate_data

    def _spatial_interpolation(self, data: ClimateData) -> ClimateData:
        logging.info("Starting spatial interpolation step...")
        # here we decide what type of spatial interpolation we need to use
        # TODO: Or, do we just want to wrap these in try/catch and see which
        #       one works? Is there a way to automatically determine which
        #       type of interpolation we need?
        try:
            logging.info("Checking spatial grid interpolation...")
            return self.interpolator.interpolate_spatial_grid(data, self.grid)
        except Exception as e:
            logging.info("Spatial grid interpolation failed.")
            logging.exception(e)

        try:
            logging.info("Checking spatial 2d interpolation...")
            return self.interpolator.interpolate_spatial_2d(data, self.grid)
        except Exception as e:
            logging.info("Spatial 2d interpolation failed.")
            logging.exception(e)

        # TODO: Is this really a warning?
        #       Some data we don't want to interpolate
        logging.warn("Returning with no spatial interpolation.")
        return data

    def _serialize(self, data: ClimateData, output_directory: str):
        logging.info("Starting serialization step...")
        # create descriptive name
        filename = data.filename + ".pkl"
        filepath = os.path.join(output_directory, filename)
        logging.info("Serializing data to %s", filepath)

        data.serialize(filepath)

    def plot_data_on_map(self, data: ClimateData, output_directory: str,
                         timestep: int = 0):
        # TODO: Plotting just the first timestep for now
        # Maybe we want to plot a video of all timesteps instead
        plt.imshow(
            data.data[timestep],
            extent=[
                data.longitudes.min(),
                data.longitudes.max(),
                data.latitudes.min(),
                data.latitudes.max()])
        plt.colorbar(fraction=0.018, pad=0.04, label="Concentration")
        formatted_date = data.dates[timestep].strftime("%Y-%m-%d")
        title = f"{data.description} {formatted_date}"
        plt.title(title)
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        filename = f"{data.description}_{formatted_date}.png"
        filepath = os.path.join(output_directory, filename)
        plt.savefig(filepath)
        plt.close()

    def preprocess(self, input_directory: str, output_directory: str,
                   only_plot_first: bool = False):
        input_filepaths = self._get_all_files_in_directory(input_directory)

        loading_strategies =\
            self._determine_loading_strategies(input_filepaths)
        if len(loading_strategies) == 0:
            logging.error("No loading strategies could be used. Exiting.")
            return

        for loading_strategy in loading_strategies:
            # load all the data for this strategy
            related_data: List[ClimateData] = []
            for i, input_file in enumerate(input_filepaths):
                # first convert the input file to the reference grid
                # NOTE: We explicitly do NOT do this for forcing variables
                logging.info("(%d/%d) Loading data for: %s",
                             i, len(input_filepaths),
                             loading_strategy.VARIABLE_KEY)
                logging.info("Input filepath: %s", input_file)
                if isinstance(loading_strategy, ForcingStrategy):
                    # load data directly from disk
                    logging.info("Not converting forcing term grid!")
                    loaded_data = loading_strategy.load_from_disk(input_file)
                else:
                    # load dataset into memory for parsing
                    variable_name = loading_strategy.VARIABLE_KEY
                    num_grid =\
                        self.cdo_wrapper.\
                        determine_grid_number_containing_variable(
                            input_file, variable_name)
                    if isinstance(loading_strategy, SicStrategy):
                        logging.info("Explicitly setting grid to reference "
                                     "before conversion")
                        dataset =\
                            self.cdo_wrapper.\
                            convert_file_and_set_grid_in_memory(
                                input_file, self.grid_filepath)
                    else:
                        dataset = self.cdo_wrapper.convert_file_grid_in_memory(
                            input_file, self.grid_filepath, num_grid)

                    # now extract data from the dataset and
                    loaded_data = loading_strategy.load_from_memory(dataset)
                    dataset.close()
                related_data.append(loaded_data)

            data = self._sort(related_data)
            data = self._temporal_interpolation(related_data)
            combined_data = self._combine(data)
            combined_data = self._spatial_interpolation(combined_data)
            logging.info("Combined data is: %s", combined_data.description)

            # put in variable specific output directory
            # skip if directory exists
            # TODO: Ideally we do this before the processing
            var_output_directory = os.path.join(
                output_directory, combined_data.description)
            if os.path.exists(var_output_directory):
                logging.warning("Skipping variable since directory exists %s",
                                var_output_directory)
                continue
            os.makedirs(var_output_directory)

            self._serialize(combined_data, var_output_directory)
            if only_plot_first:
                self.plot_data_on_map(combined_data, var_output_directory)
            else:
                for idate in range(len(combined_data.data)):
                    self.plot_data_on_map(combined_data, var_output_directory, idate)
