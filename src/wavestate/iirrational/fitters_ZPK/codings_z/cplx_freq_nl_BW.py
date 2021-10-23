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

# import scipy.linalg


class CodingFnlBW(CodingTypeZ):
    N_parameters = 2
    p_F_Hz = 0
    p_nl_BW = 0
    restabilize = False
    unstable = False

    def update(self, F_Hz, nl_BW):
        self.p_F_Hz = F_Hz
        self.p_nl_BW = nl_BW

    def reduce(self):
        return [
            self.p_F_Hz,
            self.p_nl_BW,
        ]

    def update_roots(self, r1):
        """ """
        amp = abs(r1)
        if self.restabilize:
            if amp > 1:
                print("stabilizing: ", r1)
                amp = 1 / amp
                self.lock_nl_BW = True
        self.p_F_Hz = np.angle(r1) / np.pi * self.sys.F_nyquist_Hz

        BW_Hz = (amp - 1) * self.sys.F_nyquist_Hz
        # TODO, use deadzoning instead of this method
        if BW_Hz <= 0:
            self.unstable = False
            self.p_nl_BW = (-BW_Hz) ** 0.5
        else:
            self.unstable = True
            self.p_nl_BW = BW_Hz ** 0.5
        return

    @property
    def F_Hz(self):
        return self.p_F_Hz

    @property
    def BW(self):
        if not self.unstable:
            BW = -self.p_nl_BW ** 2
        else:
            BW = self.p_nl_BW ** 2
        return BW

    def transfer(self):
        # frequency, linear nl_BW
        amp = 1 + self.BW / self.sys.F_nyquist_Hz

        r = amp * np.cos(np.pi * self.p_F_Hz / self.sys.F_nyquist_Hz)
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return (amp * amp) * Xnsq - 2 * Xn * r + 1

    def derivative(self):
        if self.disable:
            return []
        # real/imaginary part of root
        if not self.unstable:
            DampDl = -(2 * self.p_nl_BW) / self.sys.F_nyquist_Hz
            amp = 1 - self.p_nl_BW / self.sys.F_nyquist_Hz
        else:
            DampDl = +(2 * self.p_nl_BW) / self.sys.F_nyquist_Hz
            amp = 1 + self.p_nl_BW / self.sys.F_nyquist_Hz

        F_nyquist_Hz = self.sys.F_nyquist_Hz
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq

        rcos = np.cos(np.pi * self.p_F_Hz / F_nyquist_Hz)
        rsin = np.sin(np.pi * self.p_F_Hz / F_nyquist_Hz)
        return [
            (2 * amp * np.pi / self.sys.F_nyquist_Hz) * (Xn * rsin),
            ((2 * amp * DampDl) * Xnsq - Xn * (2 * rcos * DampDl)),
        ]

    def roots_c(self):
        # frequency, linear nl_BW
        amp = 1 + self.BW / self.sys.F_nyquist_Hz
        return [amp * np.exp(Ipi * abs(self.p_F_Hz) / self.sys.F_nyquist_Hz)]
