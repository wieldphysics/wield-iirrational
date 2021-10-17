# -*- coding: utf-8 -*-
"""
"""


#import numpy as np

from ..codings_cmn import (
    CodingType,
    #Ipi,
    #I2pi
)
#import scipy.linalg

class CodingFBW(CodingType):
    N_parameters = 2
    p_BW_Hz       = 0
    p_F_Hz       = 0

    def update(self, real, imag):
        self.p_BW_Hz = real
        self.p_F_Hz = imag

    def reduce(self):
        return [
            self.p_BW_Hz,
            self.p_F_Hz,
        ]

    def update_roots(self, r1):
        """
        r2, may be unspecified, in which case it is assumed to be nothing, if r1 is real, or otherwise the conjugate of r1
        """
        self.p_BW_Hz = -r1.real
        self.p_F_Hz = r1.imag
        return

    @property
    def gain_effect(self):
        if self.p_BW_Hz > 0:
            return 1
        else:
            return 1

    def transfer(self):
        #real/imaginary part of root
        r, i = -self.p_BW_Hz, self.p_F_Hz
        X    = self.sys.Xsf_grid
        Xsq  = self.sys.Xsf_grid_sq
        return ((r*r + i*i) - 2 * X * r + Xsq)

    def derivative(self):
        if self.disable:
            return []
        #real/imaginary part of root
        r, i = -self.p_BW_Hz, self.p_F_Hz
        X    = self.sys.Xsf_grid
        return [
            -2 * (r - X),
            (2 * i),
        ]

    def roots_c(self):
        #real/imaginary part of root
        r, i = -self.p_BW_Hz, abs(self.p_F_Hz)
        return [r + 1j * i]


