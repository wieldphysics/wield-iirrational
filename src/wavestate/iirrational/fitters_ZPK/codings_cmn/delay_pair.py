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


class CodingDelayPairNl(CodingType):
    N_parameters = 2
    unstable = False
    p_nlF_Hz = 0
    p_nl_BW = 0
    min_BW_Hz = 0

    def setup(self):
        self.N_parameters = 2

    def update(self, A, B):
        nlF_Hz, nl_BW = A, B
        self.p_nlF_Hz = nlF_Hz
        self.p_nl_BW = nl_BW

    def reduce(self):
        return [
            self.p_nlF_Hz,
            self.p_nl_BW,
        ]

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
    def F_Hz(self):
        return self.p_nlF_Hz ** 2

    @F_Hz.setter
    def F_Hz(self, val):
        self.p_nlF_Hz = val ** 0.5

    @property
    def BW_Hz(self):
        return self.min_BW_Hz + self.p_nl_BW ** 2

    @BW_Hz.setter
    def BW_Hz(self, val):
        if val > self.min_BW_Hz:
            self.p_nl_BW = (val - self.min_BW_Hz) ** 0.5
        else:
            self.p_nl_BW = self.min_BW_Hz ** 0.5

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
        xterm = (r * r + i * i) + Xsq
        return (xterm + 2 * r * X) / (xterm - 2 * r * X)

    def derivative(self):
        # real/imaginary part of root
        if not self.unstable:
            r = -self.BW_Hz
            DrDp = -2 * self.p_nl_BW
        else:
            r = self.BW_Hz
            DrDp = 2 * self.p_nl_BW

        i = self.F_Hz
        DiDp = 2 * self.p_nlF_Hz

        X = self.sys.Xsf_grid

        Xsq = self.sys.Xsf_grid_sq
        xterm = (r * r + i * i) + Xsq
        xfer = (xterm + 2 * r * X) / (xterm - 2 * r * X)

        return [
            DiDp * (2 * i) * xfer * (1 / (xterm + 2 * r * X) - 1 / (xterm - 2 * r * X)),
            DrDp
            * xfer
            * (
                (2 * r + 2 * X) / (xterm + 2 * r * X)
                - (2 * r - 2 * X) / (xterm - 2 * r * X)
            ),
        ]


class CodingDelayPair(CodingType):
    N_parameters = 2
    p_nlF_Hz = 0
    p_BW = 0

    def setup(self):
        self.N_parameters = 2

    def update(self, A, B):
        nlF_Hz, nl_BW = A, B
        self.p_nlF_Hz = nlF_Hz
        self.p_BW = nl_BW

    def reduce(self):
        return [
            self.p_nlF_Hz,
            self.p_BW,
        ]

    @property
    def F_Hz(self):
        return self.p_nlF_Hz ** 2

    @F_Hz.setter
    def F_Hz(self, val):
        self.p_nlF_Hz = val ** 0.5

    @property
    def BW_Hz(self):
        return self.p_BW

    @BW_Hz.setter
    def BW_Hz(self, val):
        self.p_BW = val

    def transfer(self):
        # frequency, logarithmic amplitude
        r = -self.BW_Hz
        i = self.F_Hz
        X = self.sys.Xsf_grid
        Xsq = self.sys.Xsf_grid_sq
        xterm = (r * r + i * i) + Xsq
        return (xterm + 2 * r * X) / (xterm - 2 * r * X)

    def derivative(self):
        # real/imaginary part of root
        r = -self.BW_Hz
        DrDp = -1

        i = self.F_Hz
        DiDp = 2 * self.p_nlF_Hz

        X = self.sys.Xsf_grid

        Xsq = self.sys.Xsf_grid_sq
        xterm = (r * r + i * i) + Xsq
        xfer = (xterm + 2 * r * X) / (xterm - 2 * r * X)

        return [
            DiDp * (2 * i) * xfer * (1 / (xterm + 2 * r * X) - 1 / (xterm - 2 * r * X)),
            DrDp
            * xfer
            * (
                (2 * r + 2 * X) / (xterm + 2 * r * X)
                - (2 * r - 2 * X) / (xterm - 2 * r * X)
            ),
        ]
