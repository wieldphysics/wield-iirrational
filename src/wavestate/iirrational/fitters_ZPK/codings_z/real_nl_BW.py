# -*- coding: utf-8 -*-
"""
"""


from ..codings_cmn import (
    CodingTypeZ,
)


class CodingRealnlBW(CodingTypeZ):
    N_parameters  = 1
    unstable      = False
    negative_real = False
    p_nl_BW       = 0
    min_BW_Hz     = 0

    def option_set(self, min_BW_Hz = None, **kwargs):
        super(CodingRealnlBW, self).option_set(**kwargs)
        if min_BW_Hz is not None:
            self.min_BW_Hz = min_BW_Hz
        #TODO, should modify p_nl_A for the new max_amplitude
        return

    def check_dist_limit(self, thresh = 1):
        if self.sys.distance_limit_auto >= thresh:
            self.min_BW_Hz = self.sys.distance_limit(0)
            return True
        return False

    def update(self, nl_BW):
        self.p_nl_BW = nl_BW

    def reduce(self):
        return [self.p_nl_BW]

    def update_roots(self, r1):
        assert(r1.imag == 0)
        amp = r1
        self.check_dist_limit(thresh = 1)

        if amp < 0:
            self.negative_real = True
            amp = -amp
        else:
            self.negative_real = False

        BW_Hz = (amp - 1) * self.sys.F_nyquist_Hz

        if BW_Hz <= 0:
            self.unstable = False
            BWp = -BW_Hz
        else:
            self.unstable = True
            BWp = BW_Hz
        if BWp > 2 * self.min_BW_Hz:
            self.p_nl_BW = (BWp - self.min_BW_Hz)**.5
            return True
        else:
            self.p_nl_BW = self.min_BW_Hz**.5
            return False

    @property
    def gain_effect(self):
        if self.unstable:
            return -1
        else:
            return 1

    @property
    def BW(self):
        if not self.unstable:
            BW = -self.p_nl_BW**2
        else:
            BW = self.p_nl_BW**2
        return BW

    @property
    def F_Hz(self):
        return 0

    def transfer(self):
        self.check_dist_limit(thresh = 2)
        if not self.negative_real:
            amp = 1 + self.BW / self.sys.F_nyquist_Hz
        else:
            amp = -(1 + self.BW / self.sys.F_nyquist_Hz)
        #real, log amplitude, positive
        return (1 - amp * self.sys.Xzn_grid)

    def derivative(self):
        self.check_dist_limit(thresh = 2)
        if self.disable:
            return []
        #real/imaginary part of root
        if not self.unstable:
            DampDl = -(2 * self.p_nl_BW) / self.sys.F_nyquist_Hz
            amp = 1 - self.p_nl_BW / self.sys.F_nyquist_Hz
        else:
            DampDl = +(2 * self.p_nl_BW) / self.sys.F_nyquist_Hz
            amp = 1 + self.p_nl_BW / self.sys.F_nyquist_Hz

        if self.negative_real:
            DampDl *= -1
            amp *= -1

        return [
            -DampDl * self.sys.Xzn_grid,
        ]

    def roots_r(self):
        if not self.negative_real:
            amp = 1 + self.BW / self.sys.F_nyquist_Hz
        else:
            amp = -(1 + self.BW / self.sys.F_nyquist_Hz)
        return [amp]


