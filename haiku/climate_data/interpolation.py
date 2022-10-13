
from datetime import datetime
from dateutil.relativedelta import relativedelta
from typing import List, Tuple

import numpy as np
import numpy.typing as npt


class Interpolator:

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
            f"end_time ({end_time}) is not before start_time ({start_time})"
        return (end_time.year - start_time.year) * 12 \
            + end_time.month - start_time.month
