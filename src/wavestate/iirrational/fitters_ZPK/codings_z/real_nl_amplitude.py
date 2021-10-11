# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

from ..codings_cmn import (
    CodingTypeZ,
)


class CodingRealnlA(CodingTypeZ):
    N_parameters  = 1
    unstable      = False
    p_nl_A        = 0
    max_amplitude = 1

    def update(self, nl_A):
        self.p_nl_A = nl_A

    def reduce(self):
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
            amp = 1 / amp
            self.p_nl_A = amp / (1 - amp**2)**.5
        return

    def option_set(self, minimum_BW_Hz = None, **kwargs):
        super(CodingRealnlA, self).option_set(**kwargs)
        if minimum_BW_Hz is not None:
            self.max_amplitude = 1 - (minimum_BW_Hz / self.sys.F_nyquist_Hz)
        #TODO, should modify p_nl_A for the new max_amplitude
        return

    @property
    def amplitude(self):
        if not self.unstable:
            amp = self.max_amplitude * self.p_nl_A / (1. + self.p_nl_A**2)**.5
        else:
            amp = (1. + self.p_nl_A**2)**.5 / (self.p_nl_A * self.max_amplitude)
        return amp

    @property
    def gain_effect(self):
        if self.unstable:
            return -1
        else:
            return 1

    @property
    def F_Hz(self):
        return 0

    def transfer(self):
        #real, log amplitude, positive
        return (1 - self.amplitude * self.sys.Xzn_grid)

    def derivative(self):
        if self.disable:
            return []
        #real/imaginary part of root
        if not self.unstable:
            sqp5 = (1. + self.p_nl_A**2)**.5
            #amp = self.p_nl_A / sqp5
            DampDl = (1 / sqp5 - self.p_nl_A**2 / (1. + self.p_nl_A**2)**1.5) * self.max_amplitude
        else:
            sqp5 = (1. + self.p_nl_A**2)**.5
            #amp = sqp5 / self.p_nl_A
            DampDl = (1/sqp5 - sqp5 / self.p_nl_A**2) / self.max_amplitude

        return [
            -DampDl * self.sys.Xzn_grid,
        ]

    def roots_r(self):
        return [self.amplitude]


