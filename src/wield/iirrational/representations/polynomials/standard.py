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

from ..root_bunch import root_constraints, RBAlgorithms

from .constraints import poly_constraints

fromroots = np.polynomial.polynomial.polyfromroots
val = np.polynomial.polynomial.polyval
vander = np.polynomial.polynomial.polyvander
RBalgo = RBAlgorithms(strict=False)


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
        rB = RBalgo.expect(rvec, constraint=root_constraints.mirror_imag)
    elif poly_constraints.eRoZ == constraint:
        rvec = roots(c)
        rB = RBalgo.expect(rvec, constraint=root_constraints.mirror_quad)
    return rB


def roots(c):
    return np.polynomial.polynomial.polyroots(c)


def coeff_canonicalization_gain(c):
    return c[-1]


def roots_lnG(c):
    lnG = np.log(c[-1])
    return np.polynomial.polynomial.roots(c), lnG


def fromroots_lnG(roots):
    roots = np.asarray(roots)
    return fromroots(roots), 0


def val_lnG(X, c, lnG=0):
    return val(X, c), lnG


def vander_lnG(X, N, lnG=0):
    return vander(X, N), lnG
