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

from datetime import datetime
from typing import List, Tuple
import numpy as np
import numpy.typing as npt

from haiku.training.dmdR4Haiku import DMD_R4


def compute_KMD(sorting: str, data: npt.NDArray) -> Tuple[
        npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    """Performs Koopman Mode Decomposition using DDMD-R4.
    Outputs are the eigenvalues Lambda and the normalized modes Vtn.
    """
    n, m = data.shape
    tol = n * np.finfo(float).eps  # error tolerance

    # DDMD-R4 method
    k_mode = -1
    Vtn, Lambda, rez, RQ_ref, RSI, \
        Z_ref, rez_ref, U2, AxU_k = \
        DMD_R4(data[:, :m-1], data[:, 1:m], tol, k_mode)

    # compute eigenfunctions Kefun
    # compute eigenvectors of the adjoint
    Vtn_adjoint = np.linalg.pinv(Vtn).T
    Kefun = np.zeros((m-1, m), dtype=np.complex128)
    for ifun in range(m-1):
        for isnap in range(m):
            # the kth eigenfunction evaluated at x, g_k(x) = <x, w_k>
            # where w_k is the kth eigenvector of the adjoint
            # and < , > denotes the complex inner product
            Kefun[ifun, isnap] = np.dot(data[:, isnap],
                                        np.conjugate(Vtn_adjoint[:, ifun]))

    # choose Koopman mode sorting method
    if sorting == 'mode_power':
        # sort modes and eigenvalues by mode power (large to small)
        mode_imp = np.linalg.norm(Vtn, axis=0)
    elif sorting == 'eigenfunction':
        # sort modes and eigenvalues by mean eigenfunction magnitude
        mode_imp = np.mean(np.absolute(Kefun), axis=1)

    # sort indices from large to small
    inds = np.argsort(mode_imp)
    inds = np.flip(inds, axis=0)
    mode_imp = mode_imp[inds]

    # sort the Koopman modes, eigenvalues, and eigenfunctions
    Vtn = Vtn[:, inds]
    Lambda = Lambda[inds]
    Kefun = Kefun[inds, :]

    return Lambda, Vtn, Kefun, mode_imp


def identify_mean_modes(reomega: npt.NDArray, imomega: npt.NDArray,
                        mode_imp: npt.NDArray) -> npt.NDArray:
    """Determine mean modes from provided Koopman data."""
    modeNums = np.arange(1, len(np.squeeze(mode_imp))+1)

    # allowed variance in eigenvalues
    maxReomegaDist = 0.001
    maxImomegaDist = 0.001

    # find all modes with corresponding eigenvalues within the allowed variance
    inds = np.argwhere((np.absolute(reomega) < maxReomegaDist) &
                       (np.absolute(imomega) < maxImomegaDist))
    modeNums_sel = modeNums[inds[:, 0]]

    return modeNums_sel


def identify_annual_variation_modes(reomega: npt.NDArray, imomega: npt.NDArray,
                                    mode_imp: npt.NDArray) -> npt.NDArray:
    """Determine annual variation modes from provided Koopman data."""
    modeNums = np.arange(1, len(np.squeeze(mode_imp))+1)

    # allowed variance in eigenvalues
    maxReomegaDist = 0.01
    maxImomegaDist = 0.001

    # find all modes with corresponding eigenvalues within the allowed variance
    inds = np.argwhere((np.absolute(reomega) < maxReomegaDist) &
                       (np.absolute(1/12 - np.absolute(imomega))
                        < maxImomegaDist))
    modeNums_sel = modeNums[inds[:, 0]]

    return modeNums_sel
