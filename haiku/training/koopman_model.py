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
class KoopmanModel:

    # koopman model data needed for prediction
    training_data: npt.NDArray
    vtn_real: npt.NDArray
    vtn_imaginary: npt.NDArray
    relambda: npt.NDArray
    imlambda: npt.NDArray
    reomega: npt.NDArray
    imomega: npt.NDArray
    modes: npt.NDArray
    mode_imp: npt.NDArray
    kefun: npt.NDArray

    # geographic information needed for plotting
    latitudes: npt.NDArray
    longitudes: npt.NDArray
    masked_latitudes: npt.NDArray
    masked_longitudes: npt.NDArray

    # information about training dataset
    data_source: str
    observable_name: str
    hemisphere: str
    data_frequency: str
    ensemble_number: int

    # training parameters used for prediction
    concentration_to_consider_zero: int
    round_concentration: bool
    north_lat_bound: int
    south_lat_bound: int
    number_sea_ice_observables: int
    training_start_date: int
    training_stop_date: int
