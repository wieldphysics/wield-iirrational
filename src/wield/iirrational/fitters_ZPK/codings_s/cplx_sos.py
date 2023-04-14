#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from ..codings_cmn import (
    CodingType,
    # Ipi,
    # I2pi
)


class CodingSOS(CodingType):
    N_parameters = 2
    p_c1 = 0
    p_c2 = 0
    single_root = False

    def update(self, c1, c2=None):
        if self.single_root:
            self.p_c1 = c1
        else:
            self.p_c1 = c1
            self.p_c2 = c2

    def reduce(self):
        if self.single_root:
            return [self.p_c1]
        else:
            return [self.p_c1, self.p_c2]

    def transfer(self):
        if self.single_root:
            # second order sections (2x roots either real or complex conj)
            X = self.sys.Xsf_grid
            return self.p_c1 + X
        else:
            # second order sections (2x roots either real or complex conj)
            X = self.sys.Xsf_grid
            Xsq = self.sys.Xsf_grid_sq
            return self.p_c2 + X * self.p_c1 + Xsq

    def derivative(self):
        if self.single_root:
            return [1]
        else:
            X = self.sys.Xsf_grid
            # Xsq = self.sys.Xsf_grid_sq
            return [X, 1]

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

    def update_roots(self, r1, r2=None):
        """
        r2, may be unspecified, in which case it is assumed to be nothing, if r1 is real, or otherwise the conjugate of r1
        """
        if r2 is None and r1.imag == 0:
            self.single_root = True
            self.N_parameters = 1
            self.p_c1 = -r1.real
        else:
            self.single_root = False
            self.N_parameters = 2
            if r1.imag != 0:
                # TODO check conjugates
                self.p_c2 = r1.real ** 2 + r1.imag ** 2
                self.p_c1 = -2 * r1.real
            else:
                self.p_c2 = r1 * r2
                self.p_c1 = -(r1 + r2)
        return

    def roots(self):
        # second order sections (2x roots either real or complex conj)
        if self.single_root:
            return [-self.p_c1]
        else:
            # a = self.p_c2, b = self.p_c1, c = 1
            disc = self.p_c1 * self.p_c1 - 4 * self.p_c2
            if disc >= 0:
                sqrt_disc = disc ** 0.5
                if self.p_c1 < 0:
                    r1 = (-self.p_c1 + sqrt_disc) / 2
                else:
                    r1 = (-self.p_c1 - sqrt_disc) / 2
                r2 = self.p_c2 / r1
                return [r1, r2]
            else:
                sqrt_disc = (-disc) ** 0.5
                r1 = (-self.p_c1 + sqrt_disc * 1j) / 2
                return [r1, r1.conjugate()]

    def roots_r(self):
        # second order sections (2x roots either real or complex conj)
        if self.single_root:
            return [-self.p_c1]
        else:
            # a = self.p_c2, b = self.p_c1, c = 1
            disc = self.p_c1 * self.p_c1 - 4 * self.p_c2
            if disc >= 0:
                sqrt_disc = disc ** 0.5
                if self.p_c1 < 0:
                    r1 = (-self.p_c1 + sqrt_disc) / 2
                else:
                    r1 = (-self.p_c1 - sqrt_disc) / 2
                r2 = self.p_c2 / r1
                return [r1, r2]
            else:
                return []

    def roots_c(self):
        # second order sections (2x roots either real or complex conj)
        if self.single_root:
            return []
        else:
            # complex form is always numerically stable
            # a = self.p_c2, b = self.p_c1, c = 1
            disc = self.p_c1 * self.p_c1 - 4 * self.p_c2
            if disc >= 0:
                sqrt_disc = disc ** 0.5
                return []
            else:
                sqrt_disc = (-disc) ** 0.5
                r1 = (-self.p_c1 + sqrt_disc * 1j) / 2
                return [r1]
