# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
#import numpy as np

from ..codings_cmn import (
    CodingType,
)


class CodingRealBW(CodingType):
    N_parameters = 1
    unstable     = False
    p_BW_Hz      = 0
    restabilize  = False

    def update(self, BW):
        self.p_BW_Hz = BW

    def reduce(self):
        return [self.p_BW_Hz]

    def update_roots(self, r1):
        assert(r1.imag == 0)
        if self.restabilize:
            if r1.real > 0:
                r1 = -r1
        self.p_BW_Hz = -r1.real
        return

    @property
    def F_Hz(self):
        return 0

    @property
    def BW(self):
        return self.p_BW_Hz

    @property
    def gain_effect(self):
        if self.p_BW_Hz < 0:
            return -1
        else:
            return 1

    def transfer(self):
        #real, linear BW
        r = -self.p_BW_Hz
        return self.sys.Xsf_grid - r

    def derivative(self):
        return [1]

    def roots_r(self):
        r = -self.p_BW_Hz
        return [r]

