# -*- coding: utf-8 -*-
"""
"""


import numpy as np

from .base import (
    CodingType,
)
#import scipy.linalg


class CodingRealnlA(CodingType):
    N_parameters   = 1
    unstable       = False
    p_nl_A         = 0
    lock_nl_A = None
    hide_nl_A = None

    def setup(
            self,
            hide_all       = None,
            lock_nl_A      = None,
            hide_nl_A      = None,
            hide_amplitude = None,
            disable        = None,
    ):
        if hide_all:
            hide_nl_A = True
        if hide_amplitude is not None:
            self.hide_nl_A = hide_amplitude

        if lock_nl_A is not None:
            self.lock_nl_A = lock_nl_A
        if hide_nl_A is not None:
            self.hide_nl_A = hide_nl_A
        if self.hide_nl_A:
            N_parameters = 0
        else:
            N_parameters = 1

        if disable is not None:
            self.disable = disable
        if self.disable:
            N_parameters = 0

        self.N_parameters = N_parameters

    def update(self, nl_A = None):
        if self.disable:
            assert(nl_A is None)
            return
        if self.hide_nl_A:
            assert(nl_A is None)
        else:
            if not self.lock_nl_A:
                self.p_nl_A = nl_A

    def reduce(self):
        if self.disable:
            return []
        if self.hide_nl_A:
            return []
        else:
            return [self.p_nl_A]

    def update_roots(self, r1):
        assert(r1.imag == 0)
        amp = r1

        if abs(amp) < 1:
            #amp = x / (1. + x**2)**.5
            #amp**2 = x**2 / (1. + x**2)
            #amp**2 (1. + x**2) = x**2
            #amp**2  = x**2 ( 1 - amp**2)
            #amp**2 / (1 - amp**2) = x**2
            self.unstable = False
            self.p_nl_A = amp / (1 - amp**2)**.5
        else:
            #amp = (1. + x**2)**.5 / x
            #amp**2 =  (1. + x**2) / x**2
            #amp**2 x**2 = x**2 + 1
            #amp**2  = x**2 (amp**2 - 1)
            #amp**2 / (amp**2 - 1) = x**2
            self.unstable = True
            self.p_nl_A = amp / (amp**2 - 1)**.5
        return
        return

    @property
    def amplitude(self):
        if not self.unstable:
            amp = self.p_nl_A / (1. + self.p_nl_A**2)**.5
        else:
            amp = (1. + self.p_nl_A**2)**.5 / self.p_nl_A
        return amp

    @property
    def F_Hz(self):
        return 0

    def transfer(self):
        #real, log amplitude, positive
        return (1 - self.amplitude * sys.Xn_grid)

    def derivative(self):
        if self.disable:
            return []
        #real/imaginary part of root
        if not self.unstable:
            sqp5 = (1. + self.p_nl_A**2)**.5
            #amp = self.p_nl_A / sqp5
            DampDl = 1 / sqp5 - self.p_nl_A**2 / (1. + self.p_nl_A**2)**1.5
        else:
            sqp5 = (1. + self.p_nl_A**2)**.5
            #amp = sqp5 / self.p_nl_A
            DampDl = 1/sqp5 - sqp5 / self.p_nl_A**2

        if not self.hide_nl_A:
            return [
                -DampDl * sys.Xn_grid if not self.lock_nl_A else 0,
            ]
        else:
            return []

    def roots_r(self):
        return [self.amplitude]


