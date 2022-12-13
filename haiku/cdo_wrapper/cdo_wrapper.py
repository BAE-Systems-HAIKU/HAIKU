"""Collection of methods for coordinate conversions using cdo.

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
from typing import List

from cdo import Cdo
from netCDF4 import Dataset

from haiku.cdo_wrapper.grid_information import GridParser, GridType


class CdoWrapper:

    def __init__(self) -> None:
        """Construct a CdoWrapper object."""
        self.cdo = Cdo()
        if logging.root.level == logging.DEBUG:
            logging.info("Enabling CDO debug logging")
            self.cdo.debug = True

        # verify the cdo executable is found
        if not self.cdo.check():
            raise RuntimeError("CDO not found")
        else:
            logging.info("CDO executable found: %s", self.cdo.getCdo())

    def get_num_grids_in_file(self, data_fp: str) -> int:
        """Determine the number of grids in a data file.

        Args:
            data_fp: Filepath to the input dataset file on disk.

        Returns:
            int: Number of valid grids in the dataset
        """

        return int(self.cdo.ngrids(input=data_fp)[0])

    def get_var_names_for_grid(self, data_fp: str, grid_num: int) -> List[str]:
        """Extract the list of variables present in a grid.

        Args:
            data_fp: Filepath to the input dataset file on disk
            grid_num: The number grid to use in the source dataset

        Returns:
            List[str]: The variable names found in the dataset grid
        """
        num_grids = self.get_num_grids_in_file(data_fp)
        if grid_num > num_grids:
            raise RuntimeError(
                f"File has {num_grids} grids but queried for grid {grid_num}")

        # NOTE: The space before the "-select" is very important and needed!
        return self.cdo.showname(
            f" -select,gridnum={grid_num}", input=data_fp)[0].split(" ")

    def determine_grid_number_containing_variable(self, data_fp: str,
                                                  variable: str) -> int:
        """Find the variable in the dataset and return the associated grid.

        Args:
            data_fp: Filepath to the input dataset file on disk
            variable: Case-sensitive variable name to search for

        Returns:
            int: The number of the grid containing the variable
        """
        num_grids = self.get_num_grids_in_file(data_fp)
        # indexing for grids starts at 1
        for i in range(1, num_grids+1):
            grid_vars = self.get_var_names_for_grid(data_fp, i)
            if variable in grid_vars:
                logging.info("%s found in grid %d", variable, i)
                return i
        raise ValueError("Variable %s not found in dataset: %s"
                         % (variable, data_fp))

    def convert_file_grid_in_memory(self, source_fp: str, reference_fp: str,
                                    grid: int = 1) -> Dataset:
        """Convert the dataset grid and return it in memory.

        Args:
            source_fp: Filepath to dataset to convert
            reference_fp: Filepath to grid file or dataset to use as
                          reference to convert to
            grid: The number grid to convert

        Returns:
            Dataset: The source dataset converted to the reference grid
        """
        return self.cdo.remapbil(
            "%s -selgrid,%d" % (reference_fp, grid),
            input=source_fp,
            returnCdf=True,
            options=" -f nc4")

    def convert_file_and_set_grid_in_memory(self,
                                            source_fp: str,
                                            reference_fp: str) -> Dataset:
        """Convert the dataset grid and return it in memory.

        The grid is explicitly set on the input dataset before the conversion.
        This is necessary for some datasets that do not include correct
        grid descriptions. For example, the NSIDC SIC v4 dataset needs
        this setting to happen.

        Args:
            source_fp: Filepath to dataset to convert
            reference_fp: Filepath to grid file or dataset to use as
                          reference to convert to
            grid: The number grid to convert

        Returns:
            Dataset: The source dataset converted to the reference grid
        """
        # NOTE: Need to explicitly set format to nc4 here
        #       otherwise will default to just nc
        #       and cdo compatability errors will show up
        return self.cdo.remapbil(
            "%s -setgrid,%s" % (reference_fp, reference_fp),
            input=source_fp,
            returnCdf=True,
            options=" -f nc4")

    def convert_file_grid_to_output_file(self, source_fp: str,
                                         reference_fp: str,
                                         output_fp: str,
                                         compression_level: int = 9,
                                         grid: int = 1) -> None:

        """Convert the dataset grid and save it to disk.

        Args:
            source_fp: Filepath to dataset to convert
            reference_fp: Filepath to grid file or dataset to use as
                          reference to convert to
            output_fp: Filepath to save the converted dataset to
            compression_level: Value between 0-9 (0 meaning no compression
                               and 9 being the most compressed)
            grid: The number grid to convert

        Returns:
            Dataset: The source dataset converted to the reference grid
        """
        options_str = f"-z zip_{compression_level}"
        return self.cdo.remapbil(
            "%s -selgrid,%d" % (reference_fp, grid),
            options=options_str,
            input=source_fp,
            output=output_fp)

    def get_variable_gridtype_from_dataset(self, dataset_filepath: str,
                                           variable: str) -> GridType:
        # first determine the correct grid for the specified variable
        grid_num = self.determine_grid_number_containing_variable(
            dataset_filepath, variable)

        # now we can get the grid description
        # for the specific grid in the dataset
        grid_lines = self.cdo.griddes(" -selgrid,%d" % grid_num,
                                      input=dataset_filepath,
                                      returnXArray=True)

        # using the parsed lines, we use our GridParser class
        # to generate the Grid object used internally
        grid_parser = GridParser()
        return grid_parser._determine_grid_type(grid_lines)
