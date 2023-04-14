#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from ..codings_cmn import (
    CodingType,
    # Ipi,
    # I2pi
)

# import scipy.linalg


class CodingnlFnlBW(CodingType):
    N_parameters = 2
    unstable = False
    p_nlF_Hz = 0
    p_nl_BW = 0
    min_BW_Hz = 0
    F_cutoff_Hz = None

    def setup(self):
        self.N_parameters = 2

    def update(self, A=None, B=None):
        nlF_Hz, nl_BW = A, B
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
                self.min_BW_Hz, D = self.sys.distance_limit(
                    self.F_Hz, with_derivative=True
                )
            else:
                self.min_BW_Hz, D = self.sys.distance_limit(F_Hz, with_derivative=True)
            return D
        return 0

    @property
    def gain_effect(self):
        if not self.unstable:
            return 1
        else:
            return 1

    def update_roots(self, r1):
        """ """
        # TODO should honor min_BW_Hz
        F_Hz = abs(r1.imag)
        self.check_dist_limit(F_Hz=F_Hz, thresh=1)

        ret = True
        self.F_cutoff_Hz = self.sys.F_cutoff_Hz
        relF_Hz = F_Hz / (self.F_cutoff_Hz)
        if relF_Hz >= 1:
            relF_Hz = 0.999
            ret = False
        elif relF_Hz <= -1:
            relF_Hz = 0.001
            ret = False
        self.p_nlF_Hz = relF_Hz / (1 - relF_Hz ** 2) ** 0.5

        if r1.real <= 0:
            self.unstable = False
            r1p = -r1.real
        else:
            self.unstable = True
            r1p = r1.real
        if r1p > 2 * self.min_BW_Hz:
            self.p_nl_BW = (r1p - self.min_BW_Hz) ** 0.5
        else:
            self.p_nl_BW = self.min_BW_Hz ** 0.5
            ret = False
        return ret

    @property
    def relF_Hz(self):
        val = self.p_nlF_Hz / (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        return val

    @property
    def F_Hz(self):
        return self.F_cutoff_Hz * self.relF_Hz

    @property
    def BW_Hz(self):
        return self.min_BW_Hz + self.p_nl_BW ** 2

    def transfer(self):
        self.check_dist_limit(thresh=2)
        # frequency, logarithmic amplitude
        if not self.unstable:
            r = -self.BW_Hz
        else:
            r = self.BW_Hz
        i = self.F_Hz
        X = self.sys.Xsf_grid
        Xsq = self.sys.Xsf_grid_sq
        return (r * r + i * i) - 2 * X * r + Xsq

    def derivative(self):
        D_mBW = self.check_dist_limit(thresh=2)
        # real/imaginary part of root
        if not self.unstable:
            r = -self.BW_Hz
            DrDp = -2 * self.p_nl_BW
            D_mBW *= -1
        else:
            r = self.BW_Hz
            DrDp = 2 * self.p_nl_BW

        Fsqp5 = (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        # relF_Hz   = (self.p_nlF_Hz / Fsqp5)
        DrelFHzDp = 1 / Fsqp5 - self.p_nlF_Hz ** 2 / (1.0 + self.p_nlF_Hz ** 2) ** 1.5

        i = self.F_Hz
        X = self.sys.Xsf_grid

        # print("R", r)
        # print('D_mBW', D_mBW, sys.distance_limit_auto)
        # print(sys.distance_limit(self.F_Hz, with_derivative = True))
        # print("A", 2 * i)
        # print("B", D_mBW * (2 * r - 2 * X))
        # assert(False)
        # TODO, make this work and test it..
        # D_mBW = 0

        # TODO, need to test with D_mBW
        return [
            DrelFHzDp * self.F_cutoff_Hz * (2 * i + D_mBW * (2 * r - 2 * X)),
            2 * DrDp * r - 2 * X * DrDp,
        ]

    def roots_c(self):
        if not self.unstable:
            r = -self.BW_Hz
        else:
            r = self.BW_Hz
        i = self.F_cutoff_Hz * abs(self.relF_Hz)
        # real/imaginary part of root
        return [r + 1j * i]
