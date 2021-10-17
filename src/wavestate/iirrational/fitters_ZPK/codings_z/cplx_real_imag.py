# -*- coding: utf-8 -*-
"""
"""


from ..codings_cmn import (
    CodingTypeZ,
)


class CodingRI(CodingTypeZ):
    N_parameters = 2
    p_real       = 0
    p_imag       = 0

    def update(self, real, imag):
        self.p_real = real
        self.p_imag = imag

    def reduce(self):
        return [
            self.p_real,
            self.p_imag,
        ]

    @property
    def gain_effect(self):
        if self.p_real**2 + self.p_imag**2 > 1:
            return -1
        return 1

    def update_roots(self, r1, r2 = None):
        """
        r2, may be unspecified, in which case it is assumed to be nothing, if r1 is real, or otherwise the conjugate of r1
        """
        #TODO, check r2
        self.p_real = r1.real
        self.p_imag = r1.imag
        return

    def transfer(self):
        #real/imaginary part of root
        r, i = self.p_real, self.p_imag
        Xn    = self.sys.Xzn_grid
        Xnsq  = self.sys.Xzn_grid_sq
        return ((r*r + i*i) * Xnsq - 2 * Xn * r + 1)

    def derivative(self):
        #real/imaginary part of root
        r, i = self.p_real, self.p_imag
        Xn    = self.sys.Xzn_grid
        Xnsq  = self.sys.Xzn_grid_sq
        return [
            2 * (r * Xnsq - Xn),
            (2 * i) * Xnsq,
        ]

    def roots_c(self):
        #real/imaginary part of root
        return [self.p_real + 1j * abs(self.p_imag)]


