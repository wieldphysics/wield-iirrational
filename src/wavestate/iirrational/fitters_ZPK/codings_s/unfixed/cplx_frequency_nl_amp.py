# -*- coding: utf-8 -*-
"""
"""


import numpy as np

from .base import (
    CodingType,
    Ipi,
    #I2pi
)
#import scipy.linalg


class CodingFnlA(CodingType):
    N_parameters   = 2
    unstable       = False
    p_F_Hz         = 0
    p_nl_BW        = 0
    lock_F_Hz      = False
    lock_nl_A     = False
    hide_F_Hz      = False
    hide_nl_A     = False

    def setup(
            self,
            hide_all       = None,
            lock_F_Hz      = None,
            hide_F_Hz      = None,
            lock_nl_A      = None,
            hide_nl_A      = None,
            hide_amplitude = None,
            disable        = None,
    ):
        if hide_all:
            hide_F_Hz = True
            hide_nl_A = True

        if hide_amplitude is not None:
            self.hide_nl_A = hide_amplitude

        N_parameters = 2
        if lock_F_Hz is not None:
            self.lock_F_Hz = lock_F_Hz
        if hide_F_Hz is not None:
            self.hide_F_Hz = hide_F_Hz
        if self.hide_F_Hz:
            N_parameters -= 1

        if lock_nl_A is not None:
            self.lock_nl_A = lock_nl_A
        if hide_nl_A is not None:
            self.hide_nl_A = hide_nl_A
        if self.hide_nl_A:
            N_parameters -= 1

        if disable is not None:
            self.disable = disable
        if self.disable:
            N_parameters = 0

        self.N_parameters = N_parameters

    def update(self, A = None, B = None):
        if self.disable:
            assert(A is None and B is None)
            return
        if self.hide_nl_A:
            assert(B is None)
            if self.hide_F_Hz:
                assert(A is None)
            else:
                F_Hz = A
                if not self.lock_F_Hz:
                    self.p_F_Hz = F_Hz
        else:
            if self.hide_F_Hz:

                nl_A = A
                if not self.lock_nl_A:
                    self.p_nl_BW = nl_A
            else:
                F_Hz, nl_A = A, B
                if not self.lock_F_Hz:
                    self.p_F_Hz = F_Hz

                if not self.lock_nl_A:
                    self.p_nl_BW = nl_A

    def reduce(self):
        if self.disable:
            return []
        if self.hide_nl_A:
            if self.hide_F_Hz:
                return []
            else:
                return [
                    self.p_F_Hz
                ]
        else:
            if self.hide_F_Hz:
                return [
                    self.p_nl_BW
                ]
            else:
                return [
                    self.p_F_Hz,
                    self.p_nl_BW,
                ]

    def update_roots(self, r1):
        """
        """
        F_Hz  = r1.imag
        BW_Hz = r1.real
        self.p_F_Hz = F_Hz

        if BW_Hz < 0:
            #amp = x / (1. + x**2)**.5
            #amp**2 = x**2 / (1. + x**2)
            #amp**2 (1. + x**2) = x**2
            #amp**2  = x**2 ( 1 - amp**2)
            #amp**2 / (1 - amp**2) = x**2
            self.unstable = False
            self.p_nl_BW = amp / (1 - amp**2)**.5
        else:
            #amp = (1. + x**2)**.5 / x
            #amp**2 =  (1. + x**2) / x**2
            #amp**2 x**2 = x**2 + 1
            #amp**2  = x**2 (amp**2 - 1)
            #amp**2 / (amp**2 - 1) = x**2
            self.unstable = True
            self.p_nl_BW = amp / (amp**2 - 1)**.5
        return

    @property
    def F_Hz(self):
        return self.p_F_Hz

    @property
    def BW_Hz(self):
        if not self.unstable:
            amp = self.p_nl_BW / (1. + self.p_nl_BW**2)**.5
        else:
            amp = (1. + self.p_nl_BW**2)**.5 / self.p_nl_BW
        return amp

    def transfer(self):
        #frequency, logarithmic BW_Hz
        BW          = self.BW_Hz
        F           = self.F_Hz
        X            = sys.Xn_grid
        Xsq          = sys.Xn_grid_sq
        return ((r*r + i*i) * Xsq - 2 * Xn * r + 1)

    def derivative(self):
        if self.disable:
            return []
        #real/imaginary part of root
        if not self.unstable:
            sqp5 = (1. + self.p_nl_BW**2)**.5
            amp = self.p_nl_BW / sqp5
            DampDl = 1 / sqp5 - self.p_nl_BW**2 / (1. + self.p_nl_BW**2)**1.5
        else:
            sqp5 = (1. + self.p_nl_BW**2)**.5
            amp = sqp5 / self.p_nl_BW
            DampDl = 1/sqp5 - sqp5 / self.p_nl_BW**2

        F_nyquist_Hz = sys.F_nyquist_Hz
        X            = sys.Xn_grid
        Xsq          = sys.Xn_grid_sq
        if not self.hide_F_Hz:
            if not self.hide_nl_A:
                rcos     = np.cos(np.pi * self.p_F_Hz / F_nyquist_Hz)
                rsin     = np.sin(np.pi * self.p_F_Hz / F_nyquist_Hz)
                return [
                    (2 * amp * np.pi / sys.F_nyquist_Hz) * (X * rsin) if not self.lock_F_Hz else 0,
                    ((2 * amp * DampDl) * Xsq - Xn * (2 * rcos * DampDl)) if not self.lock_nl_A else 0,
                ]
            else:
                rsin     = np.sin(np.pi * self.p_F_Hz / F_nyquist_Hz)
                return [
                    (2 * amp * np.pi / sys.F_nyquist_Hz) * (X * rsin) if not self.lock_F_Hz else 0,
                ]
        else:
            if not self.hide_nl_A:
                rcos     = np.cos(np.pi * self.p_F_Hz / F_nyquist_Hz)
                return [
                    ((2 * amp * DampDl) * Xsq - Xn * (2 * rcos * DampDl)) if not self.lock_nl_A else 0,
                ]
            else:
                return []

    def roots_c(self):
        #real/imaginary part of root
        return [self.BW_Hz + 1j * self.F_Hz]


