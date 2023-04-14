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


class CodingRealAmp(CodingType):
    N_parameters = 1
    unstable = False
    p_amplitude = 0
    lock_amplitude = False
    hide_amplitude = False
    restabilize = False

    def setup(
        self,
        hide_all=None,
        lock_amplitude=None,
        hide_amplitude=None,
        disable=None,
    ):
        if hide_all:
            hide_amplitude = True

        if lock_amplitude is not None:
            self.lock_amplitude = lock_amplitude
        if hide_amplitude is not None:
            self.hide_amplitude = hide_amplitude

        if self.hide_amplitude:
            N_parameters = 0
        else:
            N_parameters = 1

        if disable is not None:
            self.disable = disable
        if self.disable:
            N_parameters = 0

        self.N_parameters = N_parameters

    def update(self, amplitude=None):
        if self.disable:
            assert amplitude is None
            return
        if self.hide_amplitude:
            assert amplitude is None
        else:
            if not self.lock_amplitude:
                self.p_amplitude = amplitude

    def reduce(self):
        if self.disable:
            return []
        if self.hide_amplitude:
            return []
        else:
            return [self.p_amplitude]

    def update_roots(self, r1):
        assert r1.imag == 0
        if self.restabilize:
            if r1 > 1:
                r1 = 1 / r1
        self.p_amplitude = r1
        return

    @property
    def F_Hz(self):
        return 0

    @property
    def amplitude(self):
        return self.p_amplitude

    def transfer(self):
        # real, linear amplitude
        return 1 - self.p_amplitude * sys.Xn_grid

    def derivative(self):
        if self.disable:
            return []
        # real/imaginary part of root
        if not self.hide_amplitude:
            jac = [
                -sys.Xn_grid if not self.lock_amplitude else 0,
            ]
            return jac
        else:
            return []

    def roots_r(self):
        return [self.p_amplitude]


class CodingRealAmpUX(CodingRealAmp):
    restabilize = True
