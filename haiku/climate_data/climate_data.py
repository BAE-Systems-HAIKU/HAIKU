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

from dataclasses import dataclass

import numpy.typing as npt


@dataclass
class ClimateData:

    data: npt.NDArray
    date_int: npt.NDArray
    latitudes: npt.NDArray
    longitudes: npt.NDArray
    time_days: npt.NDArray
    description: str
