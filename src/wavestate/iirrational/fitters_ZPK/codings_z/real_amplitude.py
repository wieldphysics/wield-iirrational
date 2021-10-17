# -*- coding: utf-8 -*-
"""
"""


from ..codings_cmn import (
    CodingTypeZ,
)


class CodingRealAmp(CodingTypeZ):
    N_parameters   = 1
    unstable       = False
    p_amplitude    = 0
    restabilize    = False

    def update(self, amplitude):
        self.p_amplitude = amplitude

    def reduce(self):
        return [self.p_amplitude]

    def update_roots(self, r1):
        assert(r1.imag == 0)
        if self.restabilize:
            if r1 > 1:
                r1 = 1/r1
        self.p_amplitude = r1
        return

    @property
    def gain_effect(self):
        if self.p_amplitude > 1:
            return -1
        else:
            return 1

    @property
    def F_Hz(self):
        return 0

    @property
    def amplitude(self):
        return self.p_amplitude

    def transfer(self):
        #real, linear amplitude
        return 1 - self.p_amplitude * self.sys.Xzn_grid

    def derivative(self):
        return [-self.sys.Xzn_grid]

    def roots_r(self):
        return [self.p_amplitude]


class CodingRealAmpUX(CodingRealAmp):
    restabilize = True
