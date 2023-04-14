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


class CodingnlFA(CodingType):
    N_parameters = 2
    unstable = False
    p_nlF_Hz = 0
    p_amplitude = 0
    lock_nlF_Hz = False
    lock_amplitude = False
    hide_nlF_Hz = False
    hide_amplitude = False
    max_amplitude = 1

    def setup(
        self,
        hide_all=None,
        lock_nlF_Hz=None,
        hide_nlF_Hz=None,
        lock_amplitude=None,
        hide_amplitude=None,
        disable=None,
    ):
        if hide_all:
            hide_nlF_Hz = True
            hide_amplitude = True

        if hide_amplitude is not None:
            self.hide_amplitude = hide_amplitude

        N_parameters = 2
        if lock_nlF_Hz is not None:
            self.lock_nlF_Hz = lock_nlF_Hz
        if hide_nlF_Hz is not None:
            self.hide_nlF_Hz = hide_nlF_Hz
        if self.hide_nlF_Hz:
            N_parameters -= 1

        if lock_amplitude is not None:
            self.lock_amplitude = lock_amplitude
        if hide_amplitude is not None:
            self.hide_amplitude = hide_amplitude
        if self.hide_amplitude:
            N_parameters -= 1

        if disable is not None:
            self.disable = disable
        if self.disable:
            N_parameters = 0

        self.N_parameters = N_parameters

    def update(self, A=None, B=None):
        if self.disable:
            assert A is None and B is None
            return
        if self.hide_amplitude:
            assert B is None
            if self.hide_nlF_Hz:
                assert A is None
            else:
                nlF_Hz = A
                if not self.lock_nlF_Hz:
                    self.p_nlF_Hz = nlF_Hz
        else:
            if self.hide_nlF_Hz:

                amplitude = A
                if not self.lock_amplitude:
                    self.p_amplitude = amplitude
            else:
                nlF_Hz, amplitude = A, B
                if not self.lock_nlF_Hz:
                    self.p_nlF_Hz = nlF_Hz

                if not self.lock_amplitude:
                    self.p_amplitude = amplitude

    def reduce(self):
        if self.disable:
            return []
        if self.hide_amplitude:
            if self.hide_nlF_Hz:
                return []
            else:
                return [self.p_nlF_Hz]
        else:
            if self.hide_nlF_Hz:
                return [self.p_amplitude]
            else:
                return [
                    self.p_nlF_Hz,
                    self.p_amplitude,
                ]

    def option_set(self, minimum_BW_Hz=None, **kwargs):
        super(CodingnlFA, self).option_set(**kwargs)
        if minimum_BW_Hz is not None:
            self.max_amplitude = 1 - (minimum_BW_Hz / self.sys.F_nyquist_Hz)
        # TODO, should modify p_nl_A for the new max_amplitude
        return

    def update_roots(self, r1):
        """ """
        # TODO should honor max_amplitude
        amp = abs(r1)
        F_Hz = np.angle(r1) / np.pi * sys.F_nyquist_Hz
        relF_Hz = F_Hz / sys.F_cutoff_Hz
        if relF_Hz >= 1:
            relF_Hz = 0.999
        elif relF_Hz <= 0:
            relF_Hz = 0.001
        self.p_nlF_Hz = relF_Hz / (1 - relF_Hz ** 2) ** 0.5

        if amp < 1:
            # amp = x / (1. + x**2)**.5
            # amp**2 = x**2 / (1. + x**2)
            # amp**2 (1. + x**2) = x**2
            # amp**2  = x**2 ( 1 - amp**2)
            # amp**2 / (1 - amp**2) = x**2
            self.unstable = False
            self.p_amplitude = amp / (1 - amp ** 2) ** 0.5
        else:
            # amp = (1. + x**2)**.5 / x
            # amp**2 =  (1. + x**2) / x**2
            # amp**2 x**2 = x**2 + 1
            # amp**2  = x**2 (amp**2 - 1)
            # amp**2 / (amp**2 - 1) = x**2
            self.unstable = True
            self.p_amplitude = amp / (amp ** 2 - 1) ** 0.5
        return

    @property
    def relF_Hz(self):
        val = self.p_nlF_Hz / (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        return val

    @property
    def amplitude(self):
        if not self.unstable:
            amp = (
                self.max_amplitude
                * self.p_amplitude
                / (1.0 + self.p_amplitude ** 2) ** 0.5
            )
        else:
            amp = (
                (1.0 + self.p_amplitude ** 2) ** 0.5 / self.p_amplitude
            ) / self.max_amplitude
        return amp

    def transfer(self):
        # frequency, logarithmic amplitude
        amp = self.amplitude
        r = amp * np.cos(np.pi * sys.F_cutoff_Hz * self.relF_Hz / sys.F_nyquist_Hz)
        Xn = sys.Xzn_grid
        Xnsq = sys.Xzn_grid_sq
        return (amp * amp) * Xnsq - 2 * Xn * r + 1

    def derivative(self):
        if self.disable:
            return []
        # real/imaginary part of root
        if not self.unstable:
            sqp5 = (1.0 + self.p_amplitude ** 2) ** 0.5
            amp = self.max_amplitude * self.p_amplitude / sqp5
            DampDl = self.max_amplitude * (
                1 / sqp5 - self.p_amplitude ** 2 / (1.0 + self.p_amplitude ** 2) ** 1.5
            )
        else:
            sqp5 = (1.0 + self.p_amplitude ** 2) ** 0.5
            amp = sqp5 / self.p_amplitude / self.max_amplitude
            DampDl = 1 / sqp5 - sqp5 / self.p_amplitude ** 2 / self.max_amplitude

        Fsqp5 = (1.0 + self.p_nlF_Hz ** 2) ** 0.5
        relF_Hz = self.p_nlF_Hz / Fsqp5
        DrelFHzDp = 1 / Fsqp5 - self.p_nlF_Hz ** 2 / (1.0 + self.p_nlF_Hz ** 2) ** 1.5

        F_nyquist_Hz = sys.F_nyquist_Hz
        Xn = sys.Xzn_grid
        Xnsq = sys.Xzn_grid_sq
        if not self.hide_nlF_Hz:
            if not self.hide_amplitude:
                rcos = np.cos(np.pi * sys.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
                rsin = np.sin(np.pi * sys.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
                return [
                    DrelFHzDp
                    * (2 * amp * np.pi * sys.F_cutoff_Hz / sys.F_nyquist_Hz)
                    * (Xn * rsin)
                    if not self.lock_nlF_Hz
                    else 0,
                    ((2 * amp * DampDl) * Xnsq - Xn * (2 * rcos * DampDl))
                    if not self.lock_amplitude
                    else 0,
                ]
            else:
                rsin = np.sin(np.pi * sys.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
                return [
                    DrelFHzDp
                    * (2 * amp * np.pi * sys.F_cutoff_Hz / sys.F_nyquist_Hz)
                    * (Xn * rsin)
                    if not self.lock_nlF_Hz
                    else 0,
                ]
        else:
            if not self.hide_amplitude:
                rcos = np.cos(np.pi * sys.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
                return [
                    ((2 * amp * DampDl) * Xnsq - Xn * (2 * rcos * DampDl))
                    if not self.lock_amplitude
                    else 0,
                ]
            else:
                return []

    def roots_c(self):
        # real/imaginary part of root
        return [
            self.amplitude
            * np.exp(Ipi * abs(self.relF_Hz) * sys.F_cutoff_Hz / sys.F_nyquist_Hz)
        ]
