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


class CodingnlFnlBW(CodingType):
    N_parameters = 2
    unstable = False
    p_nlF_Hz = 0
    p_nl_BW = 0
    lock_nlF_Hz = False
    lock_nl_BW = False
    hide_nlF_Hz = False
    hide_nl_BW = False
    min_BW_Hz = 0

    def setup(
        self,
        hide_all=None,
        lock_nlF_Hz=None,
        hide_nlF_Hz=None,
        lock_nl_BW=None,
        hide_nl_BW=None,
        hide_amplitude=None,
        disable=None,
    ):
        if hide_all:
            hide_nlF_Hz = True
            hide_nl_BW = True

        if hide_amplitude is not None:
            self.hide_nl_BW = hide_amplitude

        N_parameters = 2
        if lock_nlF_Hz is not None:
            self.lock_nlF_Hz = lock_nlF_Hz
        if hide_nlF_Hz is not None:
            self.hide_nlF_Hz = hide_nlF_Hz
        if self.hide_nlF_Hz:
            N_parameters -= 1

        if lock_nl_BW is not None:
            self.lock_nl_BW = lock_nl_BW
        if hide_nl_BW is not None:
            self.hide_nl_BW = hide_nl_BW
        if self.hide_nl_BW:
            N_parameters -= 1

        if disable is not None:
            self.disable = disable
            if disable:
                N_parameters = 0

        self.N_parameters = N_parameters

    def update(self, A=None, B=None):
        if self.disable:
            assert A is None and B is None
            return
        if self.hide_nl_BW:
            assert B is None
            if self.hide_nlF_Hz:
                assert A is None
            else:
                nlF_Hz = A
                if not self.lock_nlF_Hz:
                    self.p_nlF_Hz = nlF_Hz
        else:
            if self.hide_nlF_Hz:

                nl_BW = A
                if not self.lock_nl_BW:
                    self.p_nl_BW = nl_BW
            else:
                nlF_Hz, nl_BW = A, B
                if not self.lock_nlF_Hz:
                    self.p_nlF_Hz = nlF_Hz

                if not self.lock_nl_BW:
                    self.p_nl_BW = nl_BW

    def reduce(self):
        if self.disable:
            return []
        if self.hide_nl_BW:
            if self.hide_nlF_Hz:
                return []
            else:
                return [self.p_nlF_Hz]
        else:
            if self.hide_nlF_Hz:
                return [self.p_nl_BW]
            else:
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

    def update_roots(self, r1):
        """ """
        # TODO should honor min_BW_Hz
        F_Hz = r1.imag
        relF_Hz = F_Hz / (sys.F_cutoff_Hz)
        if relF_Hz >= 1:
            relF_Hz = 0.999
        elif relF_Hz <= -1:
            relF_Hz = 0.001
        self.p_nlF_Hz = relF_Hz / (1 - relF_Hz ** 2) ** 0.5

        if r1.real <= 0:
            self.unstable = False
            self.p_nl_BW = (-r1.real) ** 0.5
        else:
            self.unstable = True
            self.p_nl_BW = (r1.real) ** 0.5
        return

    @property
    def relF_Hz(self):
        val = self.p_nlF_Hz / (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        return val

    @property
    def BW_Hz(self):
        return self.min_BW_Hz + self.p_nl_BW ** 2

    def transfer(self, sys):
        # frequency, logarithmic amplitude
        if not self.unstable:
            r = -self.BW_Hz
        else:
            r = self.BW_Hz
        i = sys.F_cutoff_Hz * self.relF_Hz
        X = sys.Xsf_grid
        Xsq = sys.Xsf_grid_sq
        return (r * r + i * i) - 2 * X * r + Xsq

    def derivative(self, sys):
        if self.disable:
            return []
        # real/imaginary part of root
        if not self.unstable:
            r = -self.BW_Hz
            DrDp = -2 * self.p_nl_BW
        else:
            r = self.BW_Hz
            DrDp = 2 * self.p_nl_BW

        Fsqp5 = (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        relF_Hz = self.p_nlF_Hz / Fsqp5
        DrelFHzDp = 1 / Fsqp5 - self.p_nlF_Hz ** 2 / (1.0 + self.p_nlF_Hz ** 2) ** 1.5
        i = sys.F_cutoff_Hz * relF_Hz

        X = sys.Xsf_grid
        if not self.hide_nlF_Hz:
            if not self.hide_nl_BW:
                return [
                    DrelFHzDp * sys.F_cutoff_Hz * 2 * i if not self.lock_nlF_Hz else 0,
                    2 * DrDp * r - 2 * X * DrDp if not self.lock_nl_BW else 0,
                ]
            else:
                return [
                    DrelFHzDp * sys.F_cutoff_Hz * 2 * i if not self.lock_nlF_Hz else 0,
                ]
        else:
            if not self.hide_nl_BW:
                return [
                    2 * DrDp * r - 2 * X * DrDp if not self.lock_nl_BW else 0,
                ]
            else:
                return []

    def roots_c(self):
        if not self.unstable:
            r = -self.BW_Hz
        else:
            r = self.BW_Hz
        i = sys.F_cutoff_Hz * self.relF_Hz
        # real/imaginary part of root
        return [r + 1j * i]
