#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Calculation of the RMS of a filter with unit white noise passing through a ZPK.
This is done using an integral and residue calculus. The filter must have more
P1 than Z, unless the ZPK is in the Z domain,
where there is a natural cutoff frequency.
"""

import numpy as np

import numpy.polynomial.chebyshev as cheby
import collections


def ZPKreciprocal(ZPK):
    Z, P, K = ZPK
    return P, Z, 1/K


def ZPKscalarprod(ZPK, val):
    Z, P, K = ZPK
    return Z, P, K * val


def ZPKprod(ZPK1, ZPK2):
    Z1, P1, K1 = ZPK1
    Z2, P2, K2 = ZPK2
    return (
        tuple(Z1) + tuple(Z2),
        tuple(P1) + tuple(P2),
        K1 * K2,
    )


def ZPKdiv(ZPK1, ZPK2):
    Z1, P1, K1 = ZPK1
    Z2, P2, K2 = ZPK2
    return (
        tuple(Z1) + tuple(P2),
        tuple(P1) + tuple(Z2),
        K1 / K2,
    )


def ZPKscalardiv(ZPK, val):
    Z, P, K = ZPK
    return Z, P, K / val


def ZPKdivscalar(val, ZPK):
    Z, P, K = ZPK
    return P, Z, val / K


def ZPKscalarsum(ZPK, val):
    return ZPKsum(ZPK, ((), (), val))


def ZPKsum(ZPK1, ZPK2, scale = None):
    Z1, P1, K1 = ZPK1
    Z2, P2, K2 = ZPK2

    if scale is None:
        #TODO, check for weak high-frequency zeros?
        #TODO, use newer numpy "initial" argument?
        scale = 1/max(
            np.amax(np.concatenate([np.abs(Z1), [0]])),
            np.amax(np.concatenate([np.abs(P1), [0]])),
            np.amax(np.concatenate([np.abs(Z2), [0]])),
            np.amax(np.concatenate([np.abs(P2), [0]])),
        )

    #apply now before we rescale the poles
    P3 = tuple(P1) + tuple(P2)

    Z1 = -1j * np.asarray(Z1) * scale
    P1 = -1j * np.asarray(P1) * scale
    #Not rescaling K1, even though we should, it would have to be inverted later
    Z2 = -1j * np.asarray(Z2) * scale
    P2 = -1j * np.asarray(P2) * scale
    #Not rescaling K2, even though we should, it would have to be inverted later

    b1 = K1 * (-1j * scale)**(len(P1) - len(Z1)) * cheby.chebfromroots(Z1)
    a1 = cheby.chebfromroots(P1)

    b2 = K2 * (-1j * scale)**(len(P2) - len(Z2)) * cheby.chebfromroots(Z2)
    a2 = cheby.chebfromroots(P2)

    b3_a = cheby.chebmul(b1, a2)
    b3_b = cheby.chebmul(b2, a1)
    b3 = cheby.chebadd(b3_a, b3_b)

    #print(b1, a2)
    #print(b2, a1)
    #print(b3)

    Z3 = 1j * cheby.chebroots(b3) / scale
    K3 = b3[-1] / (-1j * scale)**(len(P3) - len(Z3)) * 2**(len(b3) - 2)

    #TODO, CHECK THIS, needs test cases
    return tuple(Z3), tuple(P3), K3

