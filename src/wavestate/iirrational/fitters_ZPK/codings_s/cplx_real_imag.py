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

class CodingRI(CodingType):
    N_parameters = 2
    p_real       = 0
    p_imag       = 0
    lock_real    = False
    lock_imag    = False
    hide_real    = False
    hide_imag    = False

    def setup(
            self,
            hide_all  = None,
            lock_real = None,
            hide_real = None,
            lock_imag = None,
            hide_imag = None,
            disable   = None,
    ):
        if hide_all:
            hide_real = True
            hide_imag = True

        N_parameters = 2
        if lock_real is not None:
            self.lock_real = lock_real
        if hide_real is not None:
            self.hide_real = hide_real
        if self.hide_real:
            N_parameters -= 1

        if lock_imag is not None:
            self.lock_imag = lock_imag
        if hide_imag is not None:
            self.hide_imag = hide_imag
        if self.hide_imag:
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
        if self.hide_imag:
            assert(B is None)
            if self.hide_real:
                assert(A is None)
            else:
                real = A
                if not self.lock_real:
                    self.p_real = real
        else:
            if self.hide_real:

                imag = A
                if not self.lock_imag:
                    self.p_imag = imag
            else:
                real, imag = A, B
                if not self.lock_real:
                    self.p_real = real

                if not self.lock_imag:
                    self.p_imag = imag

    def reduce(self):
        if self.disable:
            return []
        if self.hide_imag:
            if self.hide_real:
                return []
            else:
                return [
                    self.p_imag
                ]
        else:
            if self.hide_real:
                return [
                    self.p_imag
                ]
            else:
                return [
                    self.p_real,
                    self.p_imag,
                ]

    def update_roots(self, r1):
        """
        r2, may be unspecified, in which case it is assumed to be nothing, if r1 is real, or otherwise the conjugate of r1
        """
        self.p_real = r1.real
        self.p_imag = r1.imag
        return

    def transfer(self):
        #real/imaginary part of root
        r, i = self.p_real, self.p_imag
        X    = self.sys.Xsf_grid
        Xsq  = self.sys.Xsf_grid_sq
        return ((r*r + i*i) - 2 * X * r + Xsq)

    def derivative(self):
        if self.disable:
            return []
        #real/imaginary part of root
        r, i = self.p_real, self.p_imag
        X    = self.sys.Xsf_grid
        if not self.hide_real:
            if not self.hide_imag:
                return [
                    2 * (r - X) if not self.lock_real else 0,
                    (2 * i) if not self.lock_imag else 0,
                ]
            else:
                return [
                    2 * (r - X) if not self.lock_real else 0,
                ]
        else:
            if not self.hide_imag:
                return [
                    (2 * i) if not self.lock_imag else 0,
                ]
            else:
                return []

    def roots_c(self):
        #real/imaginary part of root
        return [self.p_real + 1j * abs(self.p_imag)]


