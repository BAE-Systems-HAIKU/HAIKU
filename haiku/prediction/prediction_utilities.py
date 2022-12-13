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

from copy import deepcopy
import os
from typing import Any, Dict

import numpy as np
import numpy.typing as npt

def make_predictions(prediction_horizon: int, Vtn: npt.NDArray,
                     Lambda: npt.NDArray) -> npt.NDArray:
    """Make predictions with Koopman model."""
    # forward propagate the eigenvalues (increment power)
    # number of Koopman modes
    ell = Lambda.shape[0]
    # number of snapshots for reconstruction (include initial snapshot)
    m = prediction_horizon

    lambda_pow = np.ones((ell, m), dtype=np.complex_)

    for j in range(1, m):
        lambda_pow[:, j] = Lambda * lambda_pow[:, j-1]

    # create predictions: sum(Koopman mode x coefficient x eigenvalue)
    preds = Vtn @ lambda_pow

    # take only real valued predictions
    preds = preds.real

    # set any negative values to zero
    preds[preds < 0.0] = 0.0

    # set any values greater than 100 to 100
    preds[preds > 100.0] = 100.0

    return preds
