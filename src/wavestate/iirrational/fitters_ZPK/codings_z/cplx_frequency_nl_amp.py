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

from ..codings_cmn import (
    CodingTypeZ,
    Ipi,
)


class CodingFnlA(CodingTypeZ):
    N_parameters = 2
    unstable = False
    p_F_Hz = 0
    p_nl_A = 0

    def update(self, F_Hz, nl_A):
        self.p_F_Hz = F_Hz
        self.p_nl_A = nl_A

    def reduce(self):
        return [
            self.p_F_Hz,
            self.p_nl_A,
        ]

    def update_roots(self, r1):
        """ """
        amp = abs(r1)
        self.p_F_Hz = np.angle(r1) / np.pi * self.sys.F_nyquist_Hz

        if amp < 1:
            # amp = x / (1. + x**2)**.5
            # amp**2 = x**2 / (1. + x**2)
            # amp**2 (1. + x**2) = x**2
            # amp**2  = x**2 ( 1 - amp**2)
            # amp**2 / (1 - amp**2) = x**2
            self.unstable = False
            self.p_nl_A = amp / (1 - amp ** 2) ** 0.5
        else:
            # amp = (1. + x**2)**.5 / x
            # amp**2 =  (1. + x**2) / x**2
            # amp**2 x**2 = x**2 + 1
            # amp**2  = x**2 (amp**2 - 1)
            # amp**2 / (amp**2 - 1) = x**2
            self.unstable = True
            self.p_nl_A = amp / (amp ** 2 - 1) ** 0.5
        return

    @property
    def F_Hz(self):
        return self.p_F_Hz

    @property
    def amplitude(self):
        if not self.unstable:
            amp = self.p_nl_A / (1.0 + self.p_nl_A ** 2) ** 0.5
        else:
            amp = (1.0 + self.p_nl_A ** 2) ** 0.5 / self.p_nl_A
        return amp

    def transfer(self):
        # frequency, logarithmic amplitude
        amp = self.amplitude
        r = amp * np.cos(np.pi * self.p_F_Hz / self.sys.F_nyquist_Hz)
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return (amp * amp) * Xnsq - 2 * Xn * r + 1

    def derivative(self):
        if self.disable:
            return []
        # real/imaginary part of root
        if not self.unstable:
            sqp5 = (1.0 + self.p_nl_A ** 2) ** 0.5
            amp = self.p_nl_A / sqp5
            DampDl = 1 / sqp5 - self.p_nl_A ** 2 / (1.0 + self.p_nl_A ** 2) ** 1.5
        else:
            sqp5 = (1.0 + self.p_nl_A ** 2) ** 0.5
            amp = sqp5 / self.p_nl_A
            DampDl = 1 / sqp5 - sqp5 / self.p_nl_A ** 2

        F_nyquist_Hz = self.sys.F_nyquist_Hz
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        rcos = np.cos(np.pi * self.p_F_Hz / F_nyquist_Hz)
        rsin = np.sin(np.pi * self.p_F_Hz / F_nyquist_Hz)
        return [
            (2 * amp * np.pi / self.sys.F_nyquist_Hz) * (Xn * rsin)
            if not self.lock_F_Hz
            else 0,
            ((2 * amp * DampDl) * Xnsq - Xn * (2 * rcos * DampDl))
            if not self.lock_nl_A
            else 0,
        ]

    def roots_c(self):
        # real/imaginary part of root
        return [self.amplitude * np.exp(Ipi * abs(self.p_F_Hz) / self.sys.F_nyquist_Hz)]
