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
)


class CodingRealnlBW(CodingType):
    N_parameters = 1
    unstable = False
    p_nl_BW = 0
    min_BW_Hz = 0

    def update(self, nl_BW):
        self.p_nl_BW = nl_BW

    def reduce(self):
        return [self.p_nl_BW]

    def update_roots(self, r1):
        assert r1.imag == 0
        self.check_dist_limit(thresh=1)

        if r1.real <= 0:
            self.unstable = False
            r1p = -r1.real
        else:
            self.unstable = True
            r1p = r1.real

        if r1p > 2 * self.min_BW_Hz:
            self.p_nl_BW = (r1p - self.min_BW_Hz) ** 0.5
            return True
        else:
            self.p_nl_BW = self.min_BW_Hz ** 0.5
            return False

    @property
    def gain_effect(self):
        if self.unstable:
            return -1
        else:
            return 1

    @property
    def F_Hz(self):
        return 0

    @property
    def BW_Hz(self):
        return self.min_BW_Hz + self.p_nl_BW ** 2

    def check_dist_limit(self, thresh=1):
        if self.sys.distance_limit_auto >= thresh:
            self.min_BW_Hz = self.sys.distance_limit(0)
            return True
        return False

    def transfer(self):
        self.check_dist_limit(thresh=2)
        # frequency, logarithmic amplitude
        if not self.unstable:
            r = -self.BW_Hz
        else:
            r = self.BW_Hz
        X = self.sys.Xsf_grid
        return X - r

    def derivative(self):
        self.check_dist_limit(thresh=2)
        # real/imaginary part of root
        if not self.unstable:
            DrDp = 2 * self.p_nl_BW
        else:
            DrDp = -2 * self.p_nl_BW

        return [DrDp]

    def roots_r(self):
        if not self.unstable:
            r = -self.BW_Hz
        else:
            r = self.BW_Hz
        return [r]
