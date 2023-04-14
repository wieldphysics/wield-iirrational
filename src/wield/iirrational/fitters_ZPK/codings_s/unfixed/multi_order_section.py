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

from .base import (
    CodingType,
    Ipi,
    # I2pi
)

# import scipy.linalg


class CodingMOS(CodingType):
    # settable!
    N_parameters = None
    N_coeffs = 0
    p_coeffs = (1,)

    def setup(
        self,
        N_coeffs=None,
        disable=None,
    ):
        if N_coeffs is not None:
            self.N_coeffs = N_coeffs
            coeffs = list(self.p_coeffs)[: self.N_parameters + 1]
            self.p_coeffs = np.array(coeffs + [0] * (1 + self.N_coeffs - len(coeffs)))

        N_parameters = self.N_coeffs

        if disable is not None:
            self.disable = disable
        if self.disable:
            N_parameters = 0

        self.N_parameters = N_parameters
        return

    def update(self, *coeffs):
        if self.disable:
            assert not coeffs
            return
        assert len(coeffs) == self.N_parameters
        self.p_coeffs = coeffs

    def reduce(self):
        if self.disable:
            return []
        return self.p_coeffs

    def update_roots(self, *rs):
        """
        r2, may be unspecified, in which case it is assumed to be nothing, if r1 is real, or otherwise the conjugate of r1
        """
        return

    def transfer(self):
        # multi-order sections (number indicates the order)
        Xn = sys.Xn_grid
        Xn_power = Xn
        value = 1
        for cn in self.p_coeffs:
            value = value + cn * Xn_power
            Xn_power = Xn_power * Xn
        return value

    def derivative(self):
        # multi-order sections (number indicates the order)
        Xn = sys.Xn_grid
        Xn_power = Xn
        jac_list = []
        for cn in self.p_coeffs:
            jac_list.append(Xn_power)
            Xn_power = Xn_power * Xn
        return jac_list

    def derivative_wtrans(self, sys):
        # multi-order sections (number indicates the order)
        Xn = sys.Xn_grid
        Xn_power = Xn
        jac_list = []
        value = 1
        for cn in self.p_coeffs:
            value = value + cn * Xn_power
            jac_list.append(Xn_power)
            Xn_power = Xn_power * Xn
        return value, jac_list

    def roots(self):
        # multi-order sections (number indicates the order)
        return np.polynomial.polynomial.polyroots([1] + self.p_coeffs)
