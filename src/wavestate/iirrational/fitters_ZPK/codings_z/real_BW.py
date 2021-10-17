# -*- coding: utf-8 -*-
"""
"""


from ..codings_cmn import (
    CodingTypeZ,
)


class CodingRealBW(CodingTypeZ):
    N_parameters = 1
    unstable     = False
    p_BW         = 0
    restabilize  = False

    def update(self, BW):
        self.p_BW = BW

    def reduce(self):
        return [self.p_BW]

    def update_roots(self, r1):
        assert(r1.imag == 0)
        if self.restabilize:
            if r1 > 1:
                r1 = 1/r1
        self.p_BW = (1 - r1) * self.sys.F_nyquist_Hz
        return

    @property
    def F_Hz(self):
        return 0

    @property
    def BW(self):
        return self.p_BW

    def transfer(self):
        #real, linear BW
        amp = 1 - self.p_BW / self.sys.F_nyquist_Hz
        return 1 - amp * self.sys.Xzn_grid

    def derivative(self):
        return [
            self.sys.Xzn_grid / self.sys.F_nyquist_Hz,
        ]

    def roots_r(self):
        amp = 1 - self.p_BW / self.sys.F_nyquist_Hz
        return [amp]

