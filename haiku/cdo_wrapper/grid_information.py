"""Classes and methods for parsing grid info and storing it.

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

from dataclasses import dataclass
from enum import Enum
import numpy as np
import numpy.typing as npt
from typing import Dict, List, Tuple


class GridType(Enum):

    CURVILINEAR = "curvilinear"
    PROJECTION = "projection"
    LONLAT = "lonlat"


@dataclass
class Grid:
    """Container for x/y coordinates representing a map."""

    y: npt.NDArray[np.float32]
    x: npt.NDArray[np.float32]
    x_name: str
    y_name: str
    grid_type: GridType

    def __repr__(self) -> str:
        return (f"(y, x)=({self.y.shape[0]}, {self.x.shape[0]}), "
                f"grid_type={self.grid_type.value}]")


class GridParser:

    def parse_cdo_gridfile(self, gridfile: str) -> Grid:
        lines = []
        with open(gridfile, 'r') as f:
            lines = f.readlines()
        return self.parse_cdo_gridfile_data(lines)

    def parse_cdo_gridfile_data(self, grid_lines: List[str]) -> Grid:
        # turn the lines into an easy to use dictionary
        grid_mapping = self._extract_mapping_from_lines(grid_lines)

        # need to know what type of grid we are working with
        if "gridtype" not in grid_mapping:
            raise ValueError("No \"gridtype\" key found in lines")
        else:
            grid_type = GridType(grid_mapping["gridtype"])

        # parse the Grid information depending on GridType
        # TODO: Different approaches based on type of grid?
        grid = self._parse_lonlat(grid_mapping)
        grid.grid_type = grid_type
        return grid

    def get_gridtype_from_gridfile(self, gridfile: str) -> GridType:
        lines = []
        with open(gridfile, 'r') as f:
            lines = f.readlines()
        return self._determine_grid_type(lines)

    def _determine_grid_type(self, grid_lines: List[str]) -> GridType:
        for line in grid_lines:
            if "gridtype" in line:
                split_line = line.split("=")
                grid_type_str = split_line[1].strip()
                grid_type = GridType(grid_type_str)
                return grid_type
        raise ValueError("No \"gridtype\" key found in lines")

    def _parse_lonlat(self, grid_mapping: Dict[str, str]) -> Grid:
        # NOTE: Assuming these required keys exist
        #       if these keys don't exist, an exception will be thrown
        xsize = int(grid_mapping["xsize"])
        ysize = int(grid_mapping["ysize"])
        xfirst = float(grid_mapping["xfirst"])
        xinc = float(grid_mapping["xinc"])
        yfirst = float(grid_mapping["yfirst"])
        yinc = float(grid_mapping["yinc"])
        xname = grid_mapping["xname"]
        yname = grid_mapping["yname"]

        grid_longitudes = np.asarray(
            [xfirst + i*xinc for i in range(xsize)])
        grid_latitudes = np.asarray(
            [yfirst + i*yinc for i in range(ysize)])

        return Grid(grid_latitudes, grid_longitudes, xname, yname,
                    GridType.LONLAT)

    def _parse_curvilinear(self, grid_mapping: Dict[str, str]) -> Grid:
        # NOTE: Assuming these required keys exist
        #       if these keys don't exist, an exception will be thrown
        xval_strings = [v.strip() for v in grid_mapping["xvals"].split(" ")]
        yval_strings = [v.strip() for v in grid_mapping["yvals"].split(" ")]
        xvals = np.asarray([float(v) for v in xval_strings])
        yvals = np.asarray([float(v) for v in yval_strings])

        return Grid(yvals, xvals, GridType.CURVILINEAR)

    def _extract_mapping_from_lines(self, grid_lines: List[str]) \
            -> Dict[str, str]:
        # remove any lines that are comments
        grid_lines =\
            list(filter(lambda l:  "#" not in l, grid_lines))

        # some values can span multiple lines in these files
        # so the first thing we do here is create a key value
        # mapping combining these multiline values into a
        # single string
        # this makes it easier to parse the values afterward
        grid_mapping = {}

        # grab first line outside loop so we can keep the
        # multiline logic
        cur_key, cur_value = self._get_key_value_pair_from_line(grid_lines[0])
        for i in range(1, len(grid_lines)):
            line = grid_lines[i]
            if "=" in line:
                # NOTE: Don't overwrite value if it's already in there
                if cur_key not in grid_mapping:
                    grid_mapping[cur_key] = cur_value
                cur_key, cur_value =\
                    self._get_key_value_pair_from_line(grid_lines[i])
            else:
                # important to make sure newlines are converted to spaces
                # to keep the value string consistent
                cur_value += " "
                cur_value += line.strip()

        return grid_mapping

    def _get_key_value_pair_from_line(self, line: str) -> Tuple[str, str]:
        split_line = line.split("=")
        key = split_line[0].strip()
        value = split_line[1].strip()
        return key, value
