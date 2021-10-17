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


class CodingFBW(CodingTypeZ):
    N_parameters = 2
    p_F_Hz       = 0
    p_BW         = 0
    restabilize  = False

    def update(self, F_Hz, BW):
        self.p_F_Hz = F_Hz
        self.p_BW = BW

    def reduce(self):
        return [
            self.p_F_Hz,
            self.p_BW,
        ]

    def update_roots(self, r1):
        amp = abs(r1)
        if self.restabilize:
            if amp > 1:
                print("stabilizing: ", r1)
                amp = 1/amp
                self.lock_BW = True
        self.p_BW = (1 - amp) * self.sys.F_nyquist_Hz
        self.p_F_Hz = np.angle(r1) / np.pi * self.sys.F_nyquist_Hz
        return

    @property
    def F_Hz(self):
        return self.p_F_Hz

    @property
    def gain_effect(self):
        #always 1 since it is complex
        if self.p_BW >= 0:
            return 1
        else:
            return 1

    def transfer(self):
        amp = 1 - self.p_BW / self.sys.F_nyquist_Hz
        r   = amp * np.cos(np.pi * self.p_F_Hz / self.sys.F_nyquist_Hz)
        Xn   = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return (amp*amp) * Xnsq - 2 * Xn * r + 1

    def derivative(self):
        if self.disable:
            return []
        amp = 1 - self.p_BW / self.sys.F_nyquist_Hz
        F_nyquist_Hz = self.sys.F_nyquist_Hz
        Xn   = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        rcos     = np.cos(np.pi * self.p_F_Hz / F_nyquist_Hz)
        rsin     = np.sin(np.pi * self.p_F_Hz / F_nyquist_Hz)
        return [
            (2 * amp * np.pi / self.sys.F_nyquist_Hz) * (Xn * rsin),
            -((2 * amp) * Xnsq - Xn * (2 * rcos)) / self.sys.F_nyquist_Hz,
        ]

    def roots_c(self):
        #frequency, linear BW
        amp = (1 - self.p_BW / self.sys.F_nyquist_Hz)
        return [amp * np.exp(Ipi * abs(self.p_F_Hz) / self.sys.F_nyquist_Hz)]
