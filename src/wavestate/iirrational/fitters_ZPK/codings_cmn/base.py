# -*- coding: utf-8 -*-
"""
"""


import numpy as np
from ... import TFmath

Ipi  = np.pi * 1j
I2pi = np.pi * 2j


class BranchCutAmbiguity(Exception):
    pass


class EmptyCopy(object):
    pass


class CodingType(object):
    disable = False
    coding_id = None

    gain_effect = 1

    def __init__(self, sys):
        #preserved through deepcopy so that codings_ignore set can be maintained
        self.coding_id = id(self)
        self.sys = sys

    def clone(self, sys):
        new = EmptyCopy()
        d = dict(self.__dict__)
        d['sys'] = sys
        new.__dict__.update(d)
        new.__class__ = self.__class__
        return new

    def roots(self):
        rc = self.roots_c()
        return self.roots_r() + rc + [r.conjugate() for r in rc]

    def roots_r(self):
        return []

    def roots_c(self):
        return []

    def roots_Sf(self):
        return self.roots()

    def roots_r_Sf(self):
        return self.roots_r()

    def roots_c_Sf(self):
        return self.roots_c()

    def update_roots_Sf(self, *rs):
        return self.update_roots(*rs)

    def option_set(self, **kwargs):
        return

    def transfer_abs_sq(self):
        #real/imaginary part of root
        return TFmath.abs_sq(self.transfer())

    def derivative_wtrans(self):
        return self.transfer(), self.derivative()

    def derivative_abs_sq_wtrans(self):
        xfer = self.transfer()
        jac = self.derivative()
        jac_abs_sq = []
        for der in jac:
            jac_abs_sq.append(2 * (der.real * xfer.real + der.imag * xfer.imag))
        return TFmath.abs_sq(xfer), jac_abs_sq

    @property
    def derivative_deadzoned(self):
        return False



class CodingTypeZ(object):
    disable = False
    coding_id = None

    gain_effect = 1

    def __init__(self, sys):
        #preserved through deepcopy so that codings_ignore set can be maintained
        self.coding_id = id(self)
        self.sys = sys

    def clone(self, sys):
        new = EmptyCopy()
        d = dict(self.__dict__)
        d['sys'] = sys
        new.__dict__.update(d)
        new.__class__ = self.__class__
        return new

    def roots(self):
        rc = self.roots_c()
        return self.roots_r() + rc + [r.conjugate() for r in rc]

    def roots_r(self):
        return []

    def roots_c(self):
        return []

    def roots_Sf(self):
        rs = []
        for r in self.roots():
            if r.imag == 0:
                if r.real > 0:
                    r_Sf = (r.real - 1) * self.sys.F_nyquist_Hz
                else:
                    raise BranchCutAmbiguity()
            else:
                F_Hz = np.angle(r) / np.pi * self.sys.F_nyquist_Hz
                amp = abs(r)
                BW = (amp - 1) * self.sys.F_nyquist_Hz
                r_Sf = BW + 1j * F_Hz
            rs.append(r_Sf)
        return rs

    def roots_r_Sf(self):
        rs = []
        for r in self.roots_r():
            if r.real > 0:
                r_Sf = (r.real - 1) * self.sys.F_nyquist_Hz
            else:
                raise BranchCutAmbiguity()
            rs.append(r_Sf)
        return rs

    def roots_c_Sf(self):
        rs = []
        for r in self.roots_c():
            F_Hz = np.angle(r) / np.pi * self.sys.F_nyquist_Hz
            amp = abs(r)
            BW = (amp - 1) * self.sys.F_nyquist_Hz
            r_Sf = BW + 1j * F_Hz
            rs.append(r_Sf)
        return rs

    def update_roots_Sf(self, *rs):
        rZs = []
        for r in rs:
            F_Hz = r.imag
            if F_Hz > self.sys.F_nyquist_Hz:
                raise BranchCutAmbiguity()
            amp = 1 + r.real / self.sys.F_nyquist_Hz
            rZ = amp * np.exp(F_Hz / self.sys.F_nyquist_Hz * np.pi * 1j)
            rZs.append(rZ)

        return self.update_roots(*rZs)

    def option_set(self, **kwargs):
        return

    def transfer_abs_sq(self):
        #real/imaginary part of root
        return TFmath.abs_sq(self.transfer())

    def derivative_wtrans(self):
        return self.transfer(), self.derivative()

    def derivative_abs_sq_wtrans(self):
        xfer = self.transfer()
        jac = self.derivative()
        jac_abs_sq = []
        for der in jac:
            jac_abs_sq.append(2 * (der.real * xfer.real + der.imag * xfer.imag))
        return TFmath.abs_sq(xfer), jac_abs_sq

    @property
    def derivative_deadzoned(self):
        return False
