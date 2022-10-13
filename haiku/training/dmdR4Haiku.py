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

DMD_R4::
    Refined Rayleigh-Ritz procedure for Data Driven Dynamic Mode Decomposition
    with data driven Residual bounds(DD-DMD-RRRR).
    For details see
    [1] Z. Drmac, I. Mezic, R. Mohr: Data driven modal decompositions: analysis
        and enhancements,
        SIAM Journal on Scientific Computing, 40(4), A2253A2285.
"""

import numpy as np
import scipy.linalg
import cmath
import math
import scipy.io as sio

def DMD_R4(X, Y, toScale = 'Ignore', toCenter = 'Ignore', k_mode=-1, tol=-1., n_ref=0, ref_select=1, target=0,overWrite = 'Donot', file_save='NoSave'):
        """
        =====.........................................................................
        Input
        =====

        X (real/complex 2D-array) ::  Data matrix. Assumed with more rows than columns.

        Y (real/complex 2D-array) ::  Data matrix, Y = A*X for some otherwise
                                    inaccessible matrix A.

        to_scale (string)         ::  Specifies whether snapshot normalization is
                                    needed. Equilibrating the snapshots in ell_2
                                    is a good idea in this setting. It does not
                                    apply if a forgetting factor is incorporated.
                                --> If to_scale == 'Scale'
                                    then X = X*D, Y = Y*D where D is diagonal
                                    matrix such that the new X has unt columns.

        to_center (string)        ::  Optional data centering. For experiments only!
                                --> If to_center == 'Center'
                                    then X = X-(1/m)*X*ones(m,1)*ones(1,m) and
                                    Y = Y-(1/m)*Y*ones(m,1)*ones(1,m). Note that
                                    after centering, again, Y = A*X.
                                --> Set to_center = 'Ignore', except for R&D
                                    experimenting.

        k_mode (integer)          ::  Truncation mode for setting numerical rank k of
                                    X, based on the singular values sigma(i).
                                --> If k_mode == -1, k is selected as
                                    k = max( i : sigma(i) >= sigma(1)*tol. )
                                --> If k_mode == -2, k is selected as
                                    k = max( i>1 : sigma(i) > sigma(i-1)*tol
                                --> If k > 0 and k <= min(m,n)
                                    the value of k_mode is then understood
                                    as the caller's request to use k=k_mode
                                    dominant left singular vectors of X.

        tol (real, >0 )           ::  Tolerance threshold for truncating the
                                    singular values in defining numerical
                                    rank of scaled data. (Only for the values
                                    k_mode = -1, k_mode = -2.)
                                    [!] the recommended value is tol=n*eps,
                                    where eps is the round-off unit.

        nobal (string)            ::  Specifies whether to switch of the balancing
                                    in the function eig() that is used to
                                    compute the eigenvalues and eigenvectors
                                    of the Rayleigh quotient.
                                --> If nobal == 'NoBalance', the balancing is
                                    not used.
                                --> Set nobal = 'Ignore' if not sure what is
                                    this about.

        n_ref (integer)           ::  Specifies how many of the Ritz pairs will be
                                    refined. If n_ref exceds k, then n_ref is
                                    set to k (=the total number of Ritz pairs).

        ref_select (integer)      ::  Specifies how to select n_ref Ritz pairs
                                    that will be refined.
                                --> If ref_select == 1,
                                    refine n_ref pairs with smallest residuals
                                --> If ref_select == 2,
                                    refine n_ref pairs closest to a given target
                                    (see the input parameter target)
                                --> If ref_select == 3,
                                    refine n_ref pairs closest to the imaginary axis
                                --> If ref_select == 4,
                                    refine  n_ref pairs closest to the unit disc
                                --> If ref_select == 5,
                                    refine n_ref pairs closest to the real axis

        target (real/complex)     ::  Specifies a point of interest when selecting
                                    which Ritz pairs to refine (see ref_select)
                                    Only the Ritz values closest to target will
                                    be refined (and their Ritz vectors).

        overwrite (string)        ::  Specifes whether the refined Ritz vectors will
                                    overite the originl ones or returned in a
                                    separate array
                                --> If overwrite == 'OverWrite' then the selected
                                    and improved Ritz vectors overwite the
                                    original ones. Otherwise, they ar returned
                                    in a separate array.

        file_save (string)        ::  Specifies whether the selected variables
                                    shold be saved in a file.
                                --> If file_save == 'NoSave', then no data is saved.
                                --> Otherwise, selected variables are saved in
                                    the file file_save (given on input).
        ======........................................................................
        Output
        ======........................................................................
        Z      (real/complex  2D-array) :: Ritz vectors, normalized so that
                                        ||Z(:,i)||_2 = 1. If the refinement is
                                        requested with overwrite option, then some
                                        of the columns of Z are refined Ritz
                                        vectors. Their indices are then listed in
                                        the integer array RSI.

        Lambda (real/complex  1D-array) :: Ritz values

        rez    (real 1D-array)          :: 2-norms of the reziduals
                                        rez(i) = ||A*Z(:,i)-Lambda(:,i)*Z(:,i)||_2

        --> The following four arrays are void if the refiement is not requested.
            Use tildas if they are not needed.

        RQ_ref (real/complex 1D-array)  :: Rayeigh quotients with refined Ritz
                                        vectors.

        RSI (integer array)             :: If n_ref > 0, RSI contains the indices
                                        of those Ritz pairs that have been refined.

        Z_ref (real/complex 2D-array)   :: Refined Ritz vectors, in the case when
                                        they do not overwrite the orginal ones.

        rez_ref (real 1D-array)         :: The residuals of the refined Ritz pairs.

        --> The following two array arguments are returned for post-procesing such as
            e.g. refinement. Use tildas if they are not needed.

        U (real/complex 2-D array)      :: The left singular vectors of the matrix X.

        AxU_k (real/complex 2-D array)  :: The underlying operator A applied to the
                                        k leading left singular vectors of X.
        ...............................................................................
        ...............................................................................
        AIMDyn
        Coded by Zlatko Drmac, Department of Mathematics, University of Zagreb.
        drmac@math.hr
        Version 1.0 ,  January  3, 2017.
        Version 1.1 , November 11, 2017.
        ...............................................................................
        ...............................................................................

        The data matrix is n-by-m and it is assumed that the number of its rows
        is much larger than the number of its columns.
        """
        n, m = np.shape(X)

        """ Optionally, the data can be centered. This is a delicate issue. Only for R&D."""
        if ( toCenter =='Center'):
    #         print('Centering data')
            X = np.copy(X) - np.dot( np.multiply((1/m),(np.dot(X,np.ones([m,1])))) , np.ones([1,m])  )
            Y = np.copy(Y) - np.dot( np.multiply((1/m),(np.dot(Y,np.ones([m,1])))) , np.ones([1,m])  )

        save_to_file = not(file_save == 'NoSave')

        if (toScale == 'Scale'):
    #         print('Scaling data')
            """ column norms of the X-data matrix """
            D = np.linalg.norm(X, 2, axis=0)
            """ remove zero column """
            JX = np.squeeze(np.where(D > 0))
            if (np.size(JX) < m):
                X = np.copy(X[:, JX])
                D = np.copy(D[JX])
                Y = np.copy(Y[:, JX])
                m = np.size(X, 1)
            # Tu et al. (On DMD - paper) begin

            """ data samples are scaled to the unit sphere and then the SVD is computed """
            U, S, V = np.linalg.svd(np.multiply(X, 1./D), full_matrices=False)
        else:
    #         print('Not scaling data')
            U, S, V = np.linalg.svd(X, full_matrices=False)

        """ if tol < 0. set tolerance to the recomended value n*eps else use the user defined one """
        eps = np.finfo(np.float).eps
        if (tol < 0.):
            tol = n*eps

        """Determine the numerical rank k of X, based on the k_mode and tol"""
        if (k_mode == -1):  # recommended
            k = np.where(S/S[0] > tol)[0][-1]+1
        elif (k_mode == -2):  # more conservative, right now for R&D purposes only
            k = np.where(S[1:m]/S[0:m-1] > tol)[0][-1]+1
        elif ((k_mode > 0) & k_mode <= min(m, n)):  # the caller inputs the rank k
            k = k_mode
        else:
            print('Bad input -> parameter k_mode - it is set to acceptable value')
            k = min(m, n)

        if (toScale == 'Scale'):
            AxU_k = np.dot(np.multiply(Y, 1./D), np.multiply(np.conj(V[:k, :]).T, 1./S[:k]))
        else:
            AxU_k = np.dot(Y, np.multiply(np.conj(V[:k, :]).T, 1./S[:k]))
        if (n_ref > 0):
            Q, R = np.linalg.qr(np.concatenate((U[:,:k], AxU_k), axis=1), mode='reduced')
            """
            The Rayleigh quotient Sk=U_k' * A * U_k, expressed without using A
            Sk = (((U(:,1:k)' * Y)*diag(1./D))*V(:,1:k) ).*(ones(k,1)*(1./S(1:k))')
            """
            Sk = np.copy(R[0:k, k:2*k])  # it is important take the new instance of submatrix
            for i in range(k):
                if (R[i, i] < 0.):
                    Sk[i, :] = (-1.)*Sk[i, :]
            if (n_ref > k):
                n_ref = k
        else:
            Sk = np.dot(np.conj(U[:, :k]).T, AxU_k)

        """
        Ritz values
        """
        Lambda, W = np.linalg.eig(Sk)
        Z = np.dot(U[:, :k], W)   # projected modes
        AxZ = np.dot(AxU_k, W)

        """ evaluating residual of each mode """
        rez = np.zeros(k)
        for j in range(k):
            rez[j] = np.linalg.norm(AxZ[:,j]-Lambda[j]*Z[:, j], 2)

        if (save_to_file):
            np.savez(file_save, U=U, S=S, V=V, k=k, AxU_k=AxU_k, Sk=Sk, Z=Z, Lambda=Lambda, rez=rez)

        if n_ref > 0:
            """ refinement requested ; select the specified pairs """
            if (ref_select == 1):
                # nref pairs with smallest residuals
                idummy = np.argsort(rez)
                RSI = idummy[:n_ref]
            elif (ref_select == 2):
                # nref pairs closest to target (provided on input)
                idummy = np.argsort(np.abs(Lambda-target))
                RSI = idummy[:n_ref]
            elif (ref_select == 3):
                # nref pairs closest to the imaginary axis
                idummy = np.sort(abs(np.real(Lambda)))
                RSI = idummy[:n_ref]
            elif (ref_select == 4):
                # nref pairs closest to the unit disc
                idummy = np.argsort(np.abs(np.abs(Lambda)-1))
                RSI = idummy[:n_ref]
            elif (ref_select == 5):
                # nref pairs closest to the real axis
                idummy = np.argsort(np.abs(np.imag(Lambda)))
                RSI = idummy[:n_ref]
            else:
                print('XY_DMD_R4 :( -> parameter <ref_select> had an illegal value')

            VR = np.zeros([k, n_ref], dtype=complex)
            RQ_ref = np.zeros(n_ref, dtype=complex)
            rez_ref = np.zeros(n_ref)

            R1 = np.copy(R[0:k, k:2*k])
            R0 = R[0:k, 0:k]
            R2 = R[k:2*k, k:2*k]
            for ii in range(n_ref):
                i = RSI[ii]
                RL = np.concatenate((R1-Lambda[i]*R0, R2), axis=0)
                UL, SL, VVL = np.linalg.svd(RL, full_matrices=False)
                vl = np.conj(VVL[-1, :])
                VR[:, ii] = vl
                rez_ref[ii] = np.linalg.norm(np.dot(RL, vl), 2)
                RQ_ref[ii] = np.dot(VVL[-1, :], np.dot(Sk, vl))
            if (overWrite == 'overwrite'):
                Z[:, RSI] = np.dot(U[:, :r], VR)
                Lambda[RSI] = RQ_ref
                rez[RSI] = rez_ref
                Z_ref = np.empty(1)
            else:
                Z_ref = np.dot(U[:, :r], VR)
        else:
            RSI = np.empty(1)
            RQ_ref = np.empty(1)
            rez_ref = np.empty(1)
            Z_ref = np.empty(1)

        U2 = np.copy(U[:,:k])

        return Z, Lambda, rez, RQ_ref, RSI, Z_ref, rez_ref, U2, AxU_k


def SLS_Spectral_Reconstruct_W_NE(Z, Lambda, F, W):
    """
    SLS_Spectral_Reconstruct_W_NE uses the normal equations (NE) to solve the
    structured weighted (W) least squares problem
    ||( F-Z*diag(x)*(Lambda^0,Lambda^1,...,Lambda^(m-1) )*diag(W)||_F --> min,
    where Lambda and x are complex vectors, W is real vector with positive
    entries, and F and Z are real or complex matrices.
    ..............................................................................
    Input:
    Z      <n x ell>    real or complex matrix (the modes)
    Lambda <ell x 1>    real or complex vector, defines Vlm = Vandermonde marix
    F      <ell x m>    real or complex data
    W      <m x 1>      real vector of positive weights

    Output:
    x      <ell x 1>   the solution of the LS problem.
    """
    Q, R = np.linalg.qr(Z, mode='reduced')
    G = np.conjugate(Q.T) @ F
    ell, m = G.shape

    Vlm  = np.ones( (ell, m), dtype=np.complex_ )
    for j in range(1, m):
        Vlm[:,j] = Lambda * Vlm[:,j-1]

    # # run when no weights are specified
    # K = ( np.conjugate(R.T) @ R) * np.conjugate( Vlm @ np.conjugate(Vlm.T) )
    # b = ( np.conjugate(Vlm) * ( np.conjugate(R.T) @ G )) @ np.ones( (m,1) )
    K = ( np.conjugate(R.T) @ R) * np.conjugate( Vlm @ np.diag( np.squeeze(np.square(W), axis=1) ) @ np.conjugate(Vlm.T) )
    b = ( np.conjugate( Vlm @ np.diag(np.squeeze(W, axis=1)) ) * ( np.conjugate(R.T) @ G @ np.diag(np.squeeze(W, axis=1)) )) @ np.ones( (m,1) )

    x = np.squeeze( np.linalg.lstsq(K, b, rcond=None)[0] )

    # create scaled basis matrix by rescaling Z
    basis = np.zeros(Z.shape, dtype=np.complex_)
    for i in range(F.shape[0]):
        # basis[i,:] = np.expand_dims(Z[i,:], axis=1) * x
        basis[i,:] = Z[i,:] * x

    return x, basis
