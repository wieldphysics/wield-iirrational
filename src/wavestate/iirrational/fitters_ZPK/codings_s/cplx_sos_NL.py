# -*- coding: utf-8 -*-
"""
"""

import numpy as np

from ..codings_cmn import (
    CodingType,
    #Ipi,
    #I2pi
)


class CodingSOSnl(CodingType):
    N_parameters = 2
    p_nl_c1 = 0
    p_nl_c2 = 0
    deriv_deadzone = 1e-4

    #TODO
    max_BW_Hz = 1e4
    max_BW_HzSq = 1e8

    #should be a positive number
    min_BW_Hz = 0

    #True false means mirror to respective side
    unstable = False
    single_root = False

    def update(
            self,
            nl_c1,
            nl_c2 = None
    ):
        if self.single_root:
            self.p_nl_c1 = nl_c1
        else:
            assert(nl_c2 is not None)
            self.p_nl_c1 = nl_c1
            self.p_nl_c2 = nl_c2

    def reduce(self):
        if self.single_root:
            return [self.p_nl_c1]
        else:
            return [self.p_nl_c1, self.p_nl_c2]

    @property
    def F_Hz(self):
        if self.single_root:
            return 0

        c1 = self.p_nl_c1**2
        c1 = c1 / (c1 / self.max_BW_Hz + 1)
        c2 = self.p_nl_c2**2
        c2 = c2 / (c2 / self.max_BW_HzSq + 1)
        disc = c1*c1 - 4 * c2
        if disc > 0:
            return 0
        else:
            return ((-disc)**.5)/2

    @property
    def gain_effect(self):
        rs = self.roots_r()
        if len(rs) == 0:
            return 1
        elif len(rs) == 1:
            return 1 if rs[0] < 0 else -1
        else:
            if rs[0] < 0:
                return 1 if rs[1] < 0 else -1
            else:
                return 1 if rs[1] > 0 else -1

    def check_dist_limit(self, F_Hz = None, thresh = 1):
        if self.sys.distance_limit_auto >= thresh:
            self.max_BW_Hz = self.sys.max_BW_Hz
            self.max_BW_HzSq = self.max_BW_Hz**2
            if F_Hz is None:
                self.min_BW_Hz, D = self.sys.distance_limit(self.F_Hz, with_derivative = True)
                #assert(self.min_BW_Hz > 0)
            else:
                self.min_BW_Hz, D = self.sys.distance_limit(F_Hz, with_derivative = True)
                #assert(self.min_BW_Hz > 0)
            return D
        return 0

    def transfer(self):
        if self.single_root:
            self.check_dist_limit(F_Hz = 0, thresh = 2)
            #second order sections (2x roots either real or complex conj)
            c1 = self.p_nl_c1**2
            c1 = c1 / (c1 / self.max_BW_Hz + 1)
            X  = self.sys.Xsf_grid
            if not self.unstable:
                return (c1 + self.min_BW_Hz) + X
            else:
                return -(c1 + self.min_BW_Hz) + X
        else:
            #second order sections (2x roots either real or complex conj)
            c1 = self.p_nl_c1**2
            c2 = self.p_nl_c2**2
            c1 = c1 / (c1 / self.max_BW_Hz + 1)
            c2 = c2 / (c2 / self.max_BW_HzSq + 1)
            disc = c1*c1 - 4 * c2
            if disc > 0:
                self.check_dist_limit(F_Hz = 0, thresh = 2)
            else:
                F_Hz = ((-disc)**.5)/2
                self.check_dist_limit(F_Hz = F_Hz, thresh = 2)
            V  = self.min_BW_Hz
            X   = self.sys.Xsf_grid
            Xsq = self.sys.Xsf_grid_sq
            #V = 0
            if self.unstable:
                xfer = (c2 + V*(V + c1)) - (X * (c1 + 2 * V) - Xsq)
            else:
                xfer = (c2 + V*(V + c1)) + (X * (c1 + 2 * V) + Xsq)
            return xfer

    @property
    def derivative_deadzoned(self):
        if abs(self.p_nl_c1) < self.deriv_deadzone:
            return True
        if not self.single_root:
            if abs(self.p_nl_c2) < self.deriv_deadzone:
                return True
        return False

    def derivative(self):
        c1 = self.p_nl_c1**2
        D1 = (c1 / self.max_BW_Hz + 1)
        pD_c1 = 2 * self.p_nl_c1 / D1 * (1 - c1 / self.max_BW_Hz)
        if self.p_nl_c1 > 0:
            if self.p_nl_c1 < self.deriv_deadzone:
                pD_c1 = 2*self.deriv_deadzone
        else:
            if self.p_nl_c1 > -self.deriv_deadzone:
                pD_c1 = -2*self.deriv_deadzone

        if self.single_root:
            self.check_dist_limit(F_Hz = 0, thresh = 2)

            if not self.unstable:
                return [pD_c1]
            else:
                return [-pD_c1]
        else:
            c1 = c1 / D1
            c2 = self.p_nl_c2**2
            D2 = (c2 / self.max_BW_HzSq + 1)
            c2 = c2 / D2

            pD_c2 = 2 * self.p_nl_c2 / D2 * (1 - c2 / self.max_BW_HzSq)
            if self.p_nl_c2 > 0:
                if self.p_nl_c2 < self.deriv_deadzone:
                    pD_c2 = 2*self.deriv_deadzone
            else:
                if self.p_nl_c2 > -self.deriv_deadzone:
                    pD_c2 = -2*self.deriv_deadzone

            disc = c1*c1 - 4 * c2
            if disc > 0:
                V_D = self.check_dist_limit(F_Hz = 0, thresh = 2)
                V_D_c1 = 0
                V_D_c2 = 0
            else:
                F_Hz = ((-disc)**.5)/2
                V_D = self.check_dist_limit(F_Hz = F_Hz, thresh = 2)
                V_D_c1 = -pD_c1 * c1 / disc * V_D
                V_D_c2 = -pD_c2 / (2 * disc) * V_D

            #TODO, the derivative inclusion of V_D_c1,2 isn't tested well, but
            #that is partially because it doesn't matter much
            V_D_c1 = 0
            V_D_c2 = 0

            V  = self.min_BW_Hz

            X = self.sys.Xsf_grid
            if not self.unstable:
                return [
                    (pD_c1 + 2 * V_D_c1) * X + 2*V*(pD_c1 + V_D_c1),
                    (2 * X * V_D_c2) + pD_c2 + V_D_c2 * (2 * V + c1),
                ]
            else:
                return [
                    -(pD_c1 + 2 * V_D_c1) * X + 2*V*(pD_c1 + V_D_c1) + pD_c1 * V,
                    -(2 * X * V_D_c2) + pD_c2 + V_D_c2 * (2 * V + c1),
                ]

    def update_roots(self, r1, r2 = None):
        """
        r2, may be unspecified, in which case it is assumed to be nothing, if r1 is real, or otherwise the conjugate of r1
        """
        #TODO, incorporate effects of max_BW_Hz
        if r2 is None and r1.imag == 0:
            self.check_dist_limit(F_Hz = 0, thresh = 1)
            self.single_root = True
            self.N_parameters = 1
            #TODO, fix for S domain
            if self.unstable is not None:
                if r1.real > 0:
                    self.unstable = True
                    r1p = r1.real
                else:
                    self.unstable = False
                    r1p = -r1.real
                if r1p > self.min_BW_Hz:
                    self.p_nl_c1 = (r1p - self.min_BW_Hz)**.5
                    return True
                else:
                    self.p_nl_c1 = 0
                    return False
        else:
            self.single_root = False
            self.N_parameters = 2
            if r1.real > 0:
                self.unstable = True
                r1p = r1.real
            else:
                self.unstable = False
                r1p = -r1.real

            if r1.imag != 0:
                F_Hz = abs(r1.imag)
                self.check_dist_limit(F_Hz = F_Hz, thresh = 1)
                if r1p > self.min_BW_Hz:
                    ret = True
                    r1p = r1p - self.min_BW_Hz
                else:
                    ret = False
                    r1p = 0
                #TODO check conjugates
                self.p_nl_c2 = (r1p**2 + r1.imag**2)**.5
                self.p_nl_c1 = (2 * r1p)**.5
            else:
                self.check_dist_limit(0, thresh = 1)
                if r2 is None:
                    r2 = r1.conjugate()
                if r1p > self.min_BW_Hz:
                    ret = True
                    r1p = r1p - self.min_BW_Hz
                else:
                    ret = False
                    r1p = 0

                if self.unstable:
                    if r2 < 0:
                        raise RuntimeError("Can't share stable/unstable roots in this coding")
                    r2p = r2.real
                else:
                    if r2 > 0:
                        raise RuntimeError("Can't share stable/unstable roots in this coding")
                    r2p = -r2.real

                if r2p > self.min_BW_Hz:
                    r2p = r2p - self.min_BW_Hz
                else:
                    ret = False
                    r2p = 0
                self.p_nl_c2 = (r1p * r2p)**.5
                self.p_nl_c1 = (r1p + r2p)**.5
            return ret
        return

    def roots(self):
        #second order sections (2x roots either real or complex conj)
        c1 = self.p_nl_c1**2
        D1 = (c1 / self.max_BW_Hz + 1)
        c1 = c1 / D1
        if self.single_root:
            if not self.unstable:
                return [-c1 - self.min_BW_Hz]
            else:
                return [c1 + self.min_BW_Hz]
        else:
            c2 = self.p_nl_c2**2
            D2 = (c2 / self.max_BW_HzSq + 1)
            c2 = c2 / D2
            #a = c2, b = c1, c = 1
            disc = c1*c1 - 4 * c2
            if disc >= 0:
                sqrt_disc = disc**.5
                if c1 < 0:
                    r1 = (-c1 + sqrt_disc)/2
                else:
                    r1 = (-c1 - sqrt_disc)/2
                r2 = c2 / r1
                if not self.unstable:
                    return [r1 - self.min_BW_Hz, r2 - self.min_BW_Hz]
                else:
                    return [-r1 + self.min_BW_Hz, -r2 + self.min_BW_Hz]
            else:
                sqrt_disc = (-disc)**.5
                if not self.unstable:
                    r1 = (-c1 + sqrt_disc * 1j)/2 - self.min_BW_Hz
                else:
                    r1 = (+c1 + sqrt_disc * 1j)/2 + self.min_BW_Hz
                return [r1, r1.conjugate()]

    def roots_r(self):
        #second order sections (2x roots either real or complex conj)
        c1 = self.p_nl_c1**2
        D1 = (c1 / self.max_BW_Hz + 1)
        c1 = c1 / D1
        if self.single_root:
            if not self.unstable:
                return [-c1 - self.min_BW_Hz]
            else:
                return [c1 + self.min_BW_Hz]
        else:
            c2 = self.p_nl_c2**2
            D2 = (c2 / self.max_BW_HzSq + 1)
            c2 = c2 / D2
            #a = c2, b = c1, c = 1
            disc = c1*c1 - 4 * c2
            if disc >= 0:
                sqrt_disc = disc**.5
                if c1 < 0:
                    r1 = (-c1 + sqrt_disc)/2
                else:
                    r1 = (-c1 - sqrt_disc)/2
                r2 = c2 / r1
                if not self.unstable:
                    return [r1 - self.min_BW_Hz, r2 - self.min_BW_Hz]
                else:
                    return [-r1 + self.min_BW_Hz, -r2 + self.min_BW_Hz]
            else:
                return []

    def roots_c(self):
        #second order sections (2x roots either real or complex conj)
        c1 = self.p_nl_c1**2
        D1 = (c1 / self.max_BW_Hz + 1)
        c1 = c1 / D1
        if self.single_root:
            return []
        else:
            c2 = self.p_nl_c2**2
            D2 = (c2 / self.max_BW_HzSq + 1)
            c2 = c2 / D2
            #a = c2, b = c1, c = 1
            disc = c1*c1 - 4 * c2
            if disc >= 0:
                return []
            else:
                sqrt_disc = (-disc)**.5
                if not self.unstable:
                    r1 = (-c1 + sqrt_disc * 1j)/2 - self.min_BW_Hz
                else:
                    r1 = (+c1 + sqrt_disc * 1j)/2 + self.min_BW_Hz
                return [r1]
