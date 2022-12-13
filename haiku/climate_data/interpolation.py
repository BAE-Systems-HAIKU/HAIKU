
import logging
from datetime import datetime
from typing import List, Tuple

import numpy as np
import numpy.typing as npt
from dateutil.relativedelta import relativedelta
from haiku.climate_data.climate_data import ClimateData
from scipy import interpolate

from haiku.cdo_wrapper.grid_information import Grid


class Interpolator:

    def interpolate_climate_data(self, data: List[ClimateData]) -> List[ClimateData]:
        L: List[ClimateData] = []
        index = 1
        while index < len(data):
            climate_data_before = data[index - 1]
            climate_data_after = data[index]

            date_before = climate_data_before.dates[-1]
            date_after = climate_data_after.dates[0]
            elapsed_months = self.months_difference(
                date_before, date_after)
            if elapsed_months > 1:
                logging.info("%d months have elapsed between %s and %s",
                             elapsed_months, date_before, date_after)
                logging.info("Interpolating the missing months...")
                # we can interpret 0 as the starting time index
                # and elapsed_months as the ending time index
                # NOTE: This is only true since we look at data month by month
                data_before = climate_data_before.data[-1]
                data_after = climate_data_after.data[0]
                missing_data = self._interpolate_missing_months(
                    data_before,
                    data_after,
                    0,
                    elapsed_months
                )

                # create new ClimateData objects based on the new data
                for i, missing_datapoint in enumerate(missing_data):
                    # +1 here because of the zero indexing...
                    new_date = date_before + relativedelta(months=1+i)

                    # expand dimensions since this is a single timestep
                    # and we need to maintain the format (time, lat, lon)
                    # for ClimateData objects
                    # here we are adding the time dimension
                    missing_datapoint = np.expand_dims(
                        missing_datapoint, axis=0)

                    interpolated_cd = ClimateData(
                        missing_datapoint,
                        np.asarray([new_date]),
                        climate_data_before.latitudes.copy(),
                        climate_data_before.longitudes.copy(),
                        climate_data_before.description
                    )
                    L.append(interpolated_cd)

                index += elapsed_months - 1
            else:
                index += 1

        return L

    def interpolate_missing_data(self, data: npt.NDArray,
                                 date_ints: npt.NDArray) \
            -> Tuple[npt.NDArray, npt.NDArray]:
        """Detect missing timesteps and fill them in with interpolated data.

        In the event where we having missing data (gaps in time)
        we can interpolate between the data we have to fill in the gaps.

        This method currently expects months to do this interpolation.

        Args:
            data: Numpy array containing NxN arrays of data
                  First dimension is timestep
                  Second dimension is data
            date_ints: Numpy array containing ints describing time
                       Timesteps should be in units of months
                       The expected format is "%Y%m%d"
                       Example: 20001201

        Returns:
            Tuple containing:
            A new version of the data Numpy array with any missing data filled.
            A new version of the dateInt Numpy array with missing timesteps.
        """
        date_index = 1
        while date_index < len(date_ints):
            date_before_index = date_index - 1
            date_before = datetime.strptime(
                str(date_ints[date_before_index]), "%Y%m%d")
            date_after = datetime.strptime(
                str(date_ints[date_index]), "%Y%m%d")
            elapsed_months = self.months_difference(date_before, date_after)
            if elapsed_months > 1:
                print("%d months have elapsed between %s and %s" %
                      (elapsed_months, date_before, date_after))
                print("Interpolating the missing months...")
                # we can interpret 0 as the starting time index
                # and elapsed_months as the ending time index
                # NOTE: This is only true since we look at data month by month
                missing_data = self._interpolate_missing_months(
                    data[date_index - 1],
                    data[date_index],
                    0,
                    elapsed_months
                )

                indices_to_insert =\
                    [date_before_index + i for i in range(1, elapsed_months)]

                # insert the new datapoint and the new date snapshot
                # into the arrays
                for i, missing_datum in enumerate(missing_data):
                    # +1 here because of the zero indexing...
                    new_date = date_before + relativedelta(months=1+i)
                    new_date_int = int(new_date.strftime("%Y%m%d"))

                    index = indices_to_insert[i]
                    data = np.insert(data, index, missing_datum, axis=0)
                    date_ints = np.insert(
                        date_ints, index, new_date_int, axis=0)

                date_index += elapsed_months - 1

            else:
                date_index += 1

        return data, date_ints

    def _interpolate_missing_months(self, closest_data_before: npt.NDArray,
                                    closest_data_after: npt.NDArray,
                                    before_time_index: int,
                                    after_time_index: int) \
            -> List[npt.NDArray]:
        """Interpolate to infer missing data between timesteps.

        Args:
            closest_data_before: Data array before the gap
            closest_data_after: Data array after the gap
            before_time_index: Time index of the data before
            after_time_index: Time index of the data after

        Returns:
            List[npt.NDArray]: Collection of missing data in order
        """
        results: List[npt.NDArray] = []
        for current_time_index in range(before_time_index+1, after_time_index):
            # compute the data
            interpolated_data = self.interpolate(
                closest_data_before,
                closest_data_after,
                before_time_index,
                after_time_index,
                current_time_index
            )

            results.append(interpolated_data)
        return results

    def interpolate(self, y1: npt.NDArray, y2: npt.NDArray,
                    x1: int, x2: int, x: int) -> npt.NDArray:
        """Compute linear interpolation to infer missing month data.

        Interpolation formula: y = y1 + ((x - x1) / (x2 - x1)) * (y2 - y1)
            y -> unknown value
            x -> known value
            (x1, y1) -> coordinates below known x value
            (x2, y2) -> coordinates after known y value

        Args:
            y1: Closest known data before missing data
            y2: Closest known data after missing data
            x1: Timestep index of closest known data before
            x2: Timestep index of closest known data after
            x: Timestep index of missing data

        Returns:
            A Numpy array representing the interpolated data.
        """
        return y1 + ((x - x1) / (x2 - x1)) * (y2 - y1)

    def months_difference(self, start_time: datetime,
                          end_time: datetime) -> int:
        """Compute how many months have elapsed between 2 dates.

        Args:
            start_time: Starting datetime
            end_time: Ending datetime

        Returns:
            Number of months between the two dates

        Raises:
            AssertionError: If end_time is not greater than start_time
        """
        assert end_time > start_time, \
            f"end_time ({end_time}) should not be before start_time ({start_time})"
        return (end_time.year - start_time.year) * 12 \
            + end_time.month - start_time.month

    def interpolate_spatial_grid(self, data: ClimateData,
                                 target_grid: Grid) -> ClimateData:
        # arrays_are_equal =\
        #    (np.array_equal(data.latitudes, target_grid.latitudes)
        #        and np.array_equal(data.longitudes, target_grid.longitudes))
        array_shape_equal =\
            (data.latitudes.shape == target_grid.y.shape
                and data.longitudes.shape == target_grid.x.shape)

        if array_shape_equal:
            logging.info("Grid Shapes match: %s. Interpolation not needed.",
                         target_grid)
            return data
        else:
            logging.info("Different lat/lon than grid. Interpolating...")

            xx, yy = np.meshgrid(target_grid.x, target_grid.y)
            lat1D = np.reshape(data.latitudes, (data.latitudes.size))
            lon1D = np.reshape(data.longitudes, (data.longitudes.size))
            dataInterp = np.zeros(
                (
                    len(data.data),
                    len(target_grid.y),
                    len(target_grid.x)
                ), dtype=float)

            # TODO: CURRENTLY A HUGE PERFORMANCE BOTTLENECK
            # loop over first dimension of data, which is time
            for itime in range(len(data.data)):
                dataCur = np.squeeze(data.data[itime, :, :])
                dataCur1D = np.reshape(dataCur, (dataCur.size))
                dataInterp[itime, :, :] = interpolate.griddata(
                    (lon1D, lat1D), dataCur1D, (xx, yy), method='nearest')

            return ClimateData(
                dataInterp,
                data.dates,
                target_grid.y,
                target_grid.x,
                data.time_days,
                data.description,
                data.date_int)

    def interpolate_spatial_2d(self, data: ClimateData,
                               target_grid: Grid) -> ClimateData:
        arrays_are_equal =\
            (np.array_equal(data.latitudes, target_grid.y)
                and np.array_equal(data.longitudes, target_grid.x))

        if arrays_are_equal:
            logging.info("Interpolation not needed.")
            return data
        else:
            logging.info("Different lat/lon than grid. Interpolating...")
            dataInterp =\
                np.zeros(
                    (
                        len(data.data),
                        len(target_grid.y),
                        len(target_grid.x)),
                    dtype=float)
            # loop over first dimension of data, which is time
            for itime in range(len(data.data)):
                f = interpolate.interp2d(
                    data.longitudes,
                    data.latitudes,
                    np.squeeze(data.data[itime, :, :]))

                dataInterp[itime, :, :] =\
                    f(target_grid.x, target_grid.y)

            return ClimateData(
                dataInterp,
                data.dates,
                target_grid.y,
                target_grid.x,
                data.time_days,
                data.description,
                data.date_int
            )
