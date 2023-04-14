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


class CodingnlFnlBW(CodingTypeZ):
    N_parameters = 2
    unstable = False
    p_nlF_Hz = 0
    p_nl_BW = 0
    min_BW_Hz = 0

    def update(self, nlF_Hz, nl_BW):
        self.p_nlF_Hz = nlF_Hz
        self.p_nl_BW = nl_BW

    def reduce(self):
        return [
            self.p_nlF_Hz,
            self.p_nl_BW,
        ]

    def option_set(self, minimum_BW_Hz=None, **kwargs):
        super(CodingnlFnlBW, self).option_set(**kwargs)
        if minimum_BW_Hz is not None:
            self.min_BW_Hz = minimum_BW_Hz
        # TODO, should modify p_nl_BW for the new min_BW_Hz
        return

    def check_dist_limit(self, F_Hz=None, thresh=1):
        if self.sys.distance_limit_auto >= thresh:
            if F_Hz is None:
                self.min_BW_Hz = self.sys.distance_limit(self.F_Hz)
            else:
                self.min_BW_Hz = self.sys.distance_limit(F_Hz)
            return True
        return False

    def update_roots(self, r1):
        """ """
        # TODO should honor min_BW_Hz
        amp = abs(r1)
        F_Hz = abs(np.angle(r1) / np.pi * self.sys.F_nyquist_Hz)
        self.check_dist_limit(F_Hz=F_Hz, thresh=1)

        ret = True
        self.F_cutoff_Hz = self.sys.F_cutoff_Hz
        relF_Hz = F_Hz / self.F_cutoff_Hz
        if relF_Hz >= 1:
            relF_Hz = 0.999
            ret = False
        elif relF_Hz <= -1:
            relF_Hz = 0.001
            ret = False
        self.p_nlF_Hz = relF_Hz / (1 - relF_Hz ** 2) ** 0.5

        BW_Hz = (amp - 1) * self.sys.F_nyquist_Hz
        if BW_Hz <= 0:
            self.unstable = False
            BWp = -BW_Hz
        else:
            self.unstable = True
            BWp = BW_Hz

        if BWp > self.min_BW_Hz:
            self.p_nl_BW = (BWp - self.min_BW_Hz) ** 0.5
        else:
            self.p_nl_BW = self.min_BW_Hz ** 0.5
            ret = False
        return ret

    @property
    def gain_effect(self):
        # always 1 since it is complex
        if not self.unstable:
            return 1
        else:
            return 1

    @property
    def relF_Hz(self):
        val = self.p_nlF_Hz / (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        return val

    @property
    def F_Hz(self):
        return self.F_cutoff_Hz * self.relF_Hz

    @property
    def BW(self):
        if not self.unstable:
            BW = -self.p_nl_BW ** 2
        else:
            BW = self.p_nl_BW ** 2
        return BW

    def transfer(self):
        self.check_dist_limit(thresh=2)
        # frequency, logarithmic amplitude
        amp = 1 + self.BW / self.sys.F_nyquist_Hz

        r = amp * np.cos(np.pi * self.F_Hz / self.sys.F_nyquist_Hz)
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return (amp * amp) * Xnsq - 2 * Xn * r + 1

    def derivative(self):
        self.check_dist_limit(thresh=2)
        # real/imaginary part of root
        if not self.unstable:
            DampDl = -(2 * self.p_nl_BW) / self.sys.F_nyquist_Hz
            amp = 1 - self.p_nl_BW ** 2 / self.sys.F_nyquist_Hz
        else:
            DampDl = +(2 * self.p_nl_BW) / self.sys.F_nyquist_Hz
            amp = 1 + self.p_nl_BW ** 2 / self.sys.F_nyquist_Hz

        Fsqp5 = (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        relF_Hz = self.p_nlF_Hz / Fsqp5
        DrelFHzDp = 1 / Fsqp5 - self.p_nlF_Hz ** 2 / (1.0 + self.p_nlF_Hz ** 2) ** 1.5

        F_nyquist_Hz = self.sys.F_nyquist_Hz
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        rcos = np.cos(np.pi * self.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
        rsin = np.sin(np.pi * self.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
        return [
            DrelFHzDp
            * (2 * amp * np.pi * self.F_cutoff_Hz / self.sys.F_nyquist_Hz)
            * (Xn * rsin),
            ((2 * amp * DampDl) * Xnsq - Xn * (2 * rcos * DampDl)),
        ]

    def roots_c(self):
        amp = 1 + self.BW / self.sys.F_nyquist_Hz
        # real/imaginary part of root
        return [
            amp
            * np.exp(Ipi * abs(self.relF_Hz) * self.F_cutoff_Hz / self.sys.F_nyquist_Hz)
        ]
