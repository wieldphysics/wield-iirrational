#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np
import scipy
import scipy.signal
import scipy.linalg


def norm1DcSq(u):
    N = u.shape[-1]
    return np.dot(
        u.reshape(*u.shape[:-2], 1, N).conjugate(),
        u.reshape(*u.shape[:-2], N, 1),
    )[..., 0, 0]


def norm1DrSq(u):
    N = u.shape[-1]
    return np.dot(
        u.reshape(*u.shape[:-2], 1, N),
        u.reshape(*u.shape[:-2], N, 1),
    )[..., 0, 0]


def QR(
    mat,
    mshadow      = None,
    qmul         = [],
    qAmul        = [],
    pivoting     = False,
    method       = 'Householder',
    overwrite    = False,
    Rexact       = False,
    zero_test    = lambda x : x == 0,
    select_pivot = None,
):
    OHAUS = 0
    OGIVE = 1

    if not pivoting:
        def do_pivot(Cidx):
            return
    else:
        if select_pivot is None:
            def select_pivot(mtrx):
                Msum = np.sum(abs(mtrx)**2, axis = -2)
                if len(Msum.shape) > 1:
                    Msum = np.amax(Msum, axis = mtrx.shape[:-1])
                return np.argmax(Msum)

        if mshadow:
            pivmat = mshadow
        else:
            pivmat = mat
        pivots = list(range(mat.shape[-1]))
        def do_pivot(Cidx):
            Cidx2 = Cidx + select_pivot(pivmat[Cidx:, Cidx:])
            if Cidx2 == Cidx:
                return
            pivots[Cidx], pivots[Cidx2] = pivots[Cidx2], pivots[Cidx]
            swap_col(mat, Cidx, Cidx2)
            if mshadow is not None:
                swap_col(mshadow, Cidx, Cidx2)
            return

    method = method.lower()
    if method == 'householder':
        otype = OHAUS
    elif method == 'givens':
        otype = OGIVE
    else:
        raise RuntimeError("Unrecognized transformation mode")

    if not overwrite:
        mat = np.copy(mat)
        if mshadow is not None:
            mshadow = np.copy(mshadow)
        for idx, c in enumerate(qmul):
            qmul[idx] = np.copy(c)
        for idx, c in enumerate(qAmul):
            qAmul[idx] = np.copy(c)

    if otype == OGIVE:
        Nmin = min(mat.shape[-2], mat.shape[-1])
        for Cidx in range(0, Nmin):
            for Ridx in range(mat.shape[0]-1, Cidx, -1):
                #create a givens rotation for Q reduction on mat
                #from
                #On Computing Givens Rotations Reliably and Efficiently
                f = mat[Ridx - 1, Cidx]
                g = mat[Ridx, Cidx]
                if zero_test(g):
                    c = 1
                    cc = 1
                    s = 0
                    sc = 0
                    r = f
                elif zero_test(f):
                    c = 0
                    cc = 0
                    r = abs(g)
                    sc = g / r
                    s = sc.conjugate()
                else:
                    fa = abs(f)
                    rSQ = fa**2 + abs(g)**2
                    fsgn = f / fa
                    rr = rSQ**0.5
                    c = fa / rr
                    s = fsgn * g.conjugate() / rr
                    r = fsgn * rr
                    sc = s.conjugate()
                    cc = c.conjugate()
                M = np.array([
                    [c,   +s],
                    [-sc, cc],
                ])
                if Rexact:
                    def applyGR(mtrx):
                        mtrx[Ridx-1:Ridx+1, Cidx:] = M @ mtrx[Ridx-1:Ridx+1, Cidx:]
                else:
                    def applyGR(mtrx):
                        mtrx[Ridx-1:Ridx+1, Cidx+1:] = M @ mtrx[Ridx-1:Ridx+1, Cidx+1:]
                        mtrx[Ridx-1, Cidx] = r
                        mtrx[Ridx, Cidx] = 0
                applyGR(mat)
                #print(u)
                #print(mat[Cidx:, Cidx])
                def applyGRfull(mtrx):
                    mtrx[Ridx-1:Ridx+1, :] = M @ mtrx[Ridx-1:Ridx+1, :]
                if mshadow is not None:
                    applyGRfull(mshadow)
                for c in qAmul:
                    applyGRfull(c)
                def applyGRfullA(mtrx):
                    #mtrx[:, Ridx-1:Ridx+1] = mtrx[:, Ridx-1:Ridx+1] @ M
                    mtrx[:, Ridx-1:Ridx+1] = mtrx[:, Ridx-1:Ridx+1] @ M.conjugate().T
                for c in qmul:
                    applyGRfullA(c)
    elif otype == OHAUS:
        Nmin = min(mat.shape[-2], mat.shape[-1])
        for Cidx in range(0, Nmin):
            do_pivot(Cidx)
            #use starts as the x vector, will be modified in place
            u = np.copy(mat[Cidx:, Cidx])
            x0 = u[0]
            xNsq = norm1DcSq(u)
            xN = xNsq**0.5
            #TODO, need a better threshold test
            if zero_test(x0):
                x0 = 0
                x0N = 1
                alpha = -xN
            else:
                x0N = abs(x0)
                alpha = -(x0 / x0N) * xN
            u[0] -= alpha
            uNsq = (2*xN*(xN + (x0**2).real / x0N))
            if zero_test(uNsq):
                continue
            #uNsqtest = norm1DcSq(u)
            #import numpy.testing
            #numpy.testing.assert_almost_equal(uNsqtest, uNsq)
            tau = 2 / uNsq

            N = u.shape[0]
            uc = u.conjugate()
            if Rexact:
                def applyHR(mtrx):
                    mtrx[Cidx:, Cidx:] -= tau * np.dot(u.reshape(N, 1), np.dot(uc.reshape(1, N), mtrx[Cidx:, Cidx:]))
            else:
                def applyHR(mtrx):
                    mtrx[Cidx:, Cidx+1:] -= tau * np.dot(u.reshape(N, 1), np.dot(uc.reshape(1, N), mtrx[Cidx:, Cidx+1:]))
                    mtrx[Cidx, Cidx] = alpha
                    mtrx[Cidx+1:, Cidx] = 0
            applyHR(mat)
            #print(u)
            #print(mat[Cidx:, Cidx])
            def applyHRfull(mtrx):
                mtrx[Cidx:, :] -= tau * np.dot(u.reshape(N, 1), np.dot(uc.reshape(1, N), mtrx[Cidx:, :]))
            if mshadow is not None:
                applyHRfull(mshadow)
            for c in qAmul:
                applyHRfull(c)
            def applyHRfullA(mtrx):
                mtrx[:, Cidx:] -= tau * np.dot(np.dot(mtrx[:, Cidx:], u.reshape(N, 1)), uc.reshape(1, N))
            for c in qmul:
                applyHRfullA(c)
    else:
        raise NotImplementedError()

    ret = (mat,)

    if mshadow:
        ret = ret + (mshadow)

    if qmul:
        ret = ret + (qmul,)

    if qAmul:
        ret = ret + (qAmul,)

    if pivoting:
        ret = ret + (pivots,)

    if len(ret) == 1:
        return ret[0]
    else:
        return ret



def DSSQZ(
    A, B, C, D, E,
    tol = 1e-9,
):
    """
    Implementation of 
    COMPUTATION  OF IRREDUCIBLE  GENERALIZED STATE-SPACE REALIZATIONS ANDRAS VARGA
    using givens rotations.

    it is very slow, but numerically stable

    TODO, add pivoting,
    TODO, make it use the U-T property on E better for speed
    TODO, make it output Q and Z to apply to aux matrices, perhaps use them on C
    """
    #from icecream import ic
    #import tabulate
    Ninputs = B.shape[1]
    Nstates = A.shape[0]
    Nconstr = A.shape[1]
    Noutput = C.shape[0]

    BA, E = scipy.linalg.qr_multiply(
        E,
        np.hstack([B, A]),
        pivoting = False,
        mode = 'left'
    )

    Nmin = min(Nconstr, Nstates)
    for CidxBA in range(0, Nmin - 1):
        for RidxBA in range(Nconstr-1, CidxBA, -1):
            #create a givens rotation for Q reduction on BA
            BAv0 = BA[RidxBA - 1, CidxBA]
            BAv1 = BA[RidxBA, CidxBA]
            BAvSq = BAv0**2 + BAv1**2
            if BAvSq < tol:
                continue
            BAvAbs = BAvSq**0.5
            c = BAv1 / BAvAbs
            s = BAv0 / BAvAbs
            M = np.array([
                [s, +c],
                [-c, s]
            ])
            BA[RidxBA-1:RidxBA+1, :] = M @ BA[RidxBA-1:RidxBA+1, :]

            #TODO, use the U-T to be more efficient
            E[RidxBA-1:RidxBA+1, :] = M @ E[RidxBA-1:RidxBA+1, :]

            Cidx = RidxBA
            Ridx = RidxBA

            #row and col swap
            Ev0 = E[Ridx, Cidx-1]
            Ev1 = E[Ridx, Cidx]
            EvSq = Ev0**2 + Ev1**2
            if EvSq < tol:
                continue
            EvAbs = EvSq**0.5
            c = Ev0 / EvAbs
            s = Ev1 / EvAbs
            MT = np.array([
                [s, +c],
                [-c, s]
            ])
            BA[:, Ninputs:][:, Cidx-1:Cidx+1] = BA[:, Ninputs:][:, Cidx-1:Cidx+1] @ MT
            C[:, Cidx-1:Cidx+1] = C[:, Cidx-1:Cidx+1] @ MT
            #TODO, use the U-T to be more efficient
            E[:, Cidx-1:Cidx+1] = E[:, Cidx-1:Cidx+1] @ MT

    B = BA[:, :Ninputs]
    A = BA[:, Ninputs:]
    return A, B, C, D, E


def swap_col(m, Cidx1, Cidx2):
    temp = np.copy(m[:, Cidx2])
    m[:, Cidx2] = m[:, Cidx1]
    m[:, Cidx1] = temp

def swap_row(m, Cidx1, Cidx2):
    temp = np.copy(m[Cidx2, :])
    m[Cidx2, :] = m[Cidx1, :]
    m[Cidx1, :] = temp

