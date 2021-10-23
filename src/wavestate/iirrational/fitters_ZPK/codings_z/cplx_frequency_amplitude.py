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


class CodingFA(CodingTypeZ):
    N_parameters = 2
    p_F_Hz = 0
    p_amplitude = 0
    restabilize = False

    def update(self, F_Hz, amplitude):
        self.p_F_Hz = F_Hz
        self.p_amplitude = amplitude

    def reduce(self):
        return [
            self.p_F_Hz,
            self.p_amplitude,
        ]

    def update_roots(self, r1):
        """ """
        amp = abs(r1)
        if self.restabilize:
            if amp > 1:
                print("stabilizing: ", r1)
                amp = 1 / amp
                self.lock_amplitude = True
        self.p_amplitude = amp
        self.p_F_Hz = np.angle(r1) / np.pi * self.sys.F_nyquist_Hz
        return

    @property
    def F_Hz(self):
        return self.p_F_Hz

    @property
    def gain_effect(self):
        # always 1 since it is complex
        if self.p_amplitude <= 1:
            return 1
        else:
            return 1

    def transfer(self):
        # frequency, linear amplitude
        amp = self.p_amplitude
        r = amp * np.cos(np.pi * self.p_F_Hz / self.sys.F_nyquist_Hz)
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return (amp * amp) * Xnsq - 2 * Xn * r + 1

    def derivative(self):
        if self.disable:
            return []
        # real/imaginary part of root
        amp = self.p_amplitude
        F_nyquist_Hz = self.sys.F_nyquist_Hz
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        rcos = np.cos(np.pi * self.p_F_Hz / F_nyquist_Hz)
        rsin = np.sin(np.pi * self.p_F_Hz / F_nyquist_Hz)
        return [
            (2 * amp * np.pi / self.sys.F_nyquist_Hz) * (Xn * rsin),
            ((2 * amp) * Xnsq - Xn * (2 * rcos)),
        ]

    def roots_c(self):
        # frequency, linear amplitude
        return [
            self.p_amplitude * np.exp(Ipi * abs(self.p_F_Hz) / self.sys.F_nyquist_Hz)
        ]


class CodingFAUX(CodingFA):
    restabilize = True
