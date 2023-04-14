#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np
import warnings
from .constraints import poly_constraints
from ..root_bunch import root_constraints, RBAlgorithms

log2 = np.log(2)

RBalgo = RBAlgorithms(
    strict=False,
    lax_line_tol=5e-1,
)


def coeff_canonicalization_gain(c):
    # print('CCG: ', c[-1])
    return c[-1]


def roots(c):
    return np.polynomial.chebyshev.chebroots(c)


def roots_lnG(c):
    lnG = np.log(c[-1])
    return np.polynomial.chebyshev.chebroots(c), lnG


def roots_rB(c, constraint):
    if poly_constraints.even_real <= constraint:
        assert np.all(c[::2].imag == 0)
    if poly_constraints.odd_real <= constraint:
        assert np.all(c[1::2].imag == 0)
    if poly_constraints.odd_imag <= constraint:
        assert np.all(c[1::2].real == 0)

    if poly_constraints.no_constraint == constraint:
        rvec = roots(c)
        rB = RBalgo.expect(rvec, constraint=root_constraints.no_constraint)
    elif poly_constraints.eRoR == constraint:
        rvec = roots(c)
        rB = RBalgo.expect(rvec, constraint=root_constraints.mirror_real)
    elif poly_constraints.eRoI == constraint:
        rvec = roots(c)
        rB = RBalgo.expect(
            rvec, constraint=root_constraints.mirror_imag, allow_unknown=True
        )
        if len(rB.u) > 0:
            warnings.warn(
                "Unmirrored root in mirror_imag polynomial root finder (need to upgrade this algorithm)"
            )
            # HACK
            # clear any unknown
            rB.u = np.array([])
    elif poly_constraints.eRoZ == constraint:
        rvec = roots(c)
        rB = RBalgo.expect(rvec, constraint=root_constraints.mirror_quad)
    return rB


def fromroots_lnG(roots):
    lnG = log2 * (len(roots) - 1)
    c = np.polynomial.chebyshev.chebfromroots(roots)
    cG = coeff_canonicalization_gain(c)
    c /= cG
    lnG += np.log(cG)
    return c, lnG


def val_lnG(X, c, lnG=0):
    # the LOG2 is because the last coefficient is assumed to be scaled to one (as is done in fromroots)
    lnG_out = -(len(c) - 2) * log2
    return np.polynomial.chebyshev.chebval(X, c), lnG + lnG_out


def vander_lnG(X, N, lnG=0):
    # the LOG2 is because the last coefficient is assumed to be scaled to one (as is done in fromroots)
    lnG_out = -(N - 1) * log2
    return np.polynomial.chebyshev.chebvander(X, N), lnG + lnG_out


def companion(c):
    return np.polynomial.chebyshev.chebcompanion(c)
