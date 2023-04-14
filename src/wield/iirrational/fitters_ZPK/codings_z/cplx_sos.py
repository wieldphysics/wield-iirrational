#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from ..codings_cmn import (
    CodingTypeZ,
)


class CodingSOS(CodingTypeZ):
    N_parameters = 2
    p_c1 = 0
    p_c2 = 0
    # none means don't enforce mirroring, True false means mirror to respective side
    unstable = None

    def update(self, c1, c2):
        self.p_c1 = c1
        self.p_c2 = c2

    def reduce(self):
        return [self.p_c1, self.p_c2]

    def effective_params(self):
        if self.unstable is None:
            return self.p_c1, self.p_c2
        elif self.unstable:
            return self.p_c1, self.p_c2
        else:
            return self.p_c1, self.p_c2

    def transfer(self):
        # second order sections (2x roots either real or complex conj)
        c1, c2 = self.effective_params()
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return c2 * Xnsq + Xn * c1 + 1

    def derivative(self):
        Xn = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return [Xn, Xnsq]

    def update_roots(self, r1, r2=None):
        """
        r2, may be unspecified, in which case it is assumed to be nothing, if r1 is real, or otherwise the conjugate of r1
        """
        if self.unstable is not None:
            if abs(r1) >= 1:
                self.unstable = True
            else:
                self.unstable = False

        if r2 is None:
            if r1.imag != 0:
                self.p_c2 = r1.real ** 2 + r1.imag ** 2
                self.p_c1 = -2 * r1.real
            else:
                self.p_c2 = 0
                self.p_c1 = -r1.real
        else:
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

        c1, c2 = self.effective_params()
        # a = c2, b = c1, c = 1
        disc = c1 * c1 - 4 * c2
        if disc >= 0:
            sqrt_disc = disc ** 0.5
            if c1 < 0:
                r1 = (-c1 + sqrt_disc) / 2
            else:
                r1 = (-c1 - sqrt_disc) / 2
            r2 = c2 / r1
            return [r1, r2]
        else:
            sqrt_disc = (-disc) ** 0.5
            r1 = (-c1 + sqrt_disc * 1j) / 2
            return [r1, r1.conjugate()]

    def roots_r(self):
        # second order sections (2x roots either real or complex conj)

        c1, c2 = self.effective_params()
        # a = c2, b = c1, c = 1
        disc = c1 * c1 - 4 * c2
        if disc >= 0:
            sqrt_disc = disc ** 0.5
            if c1 < 0:
                r1 = (-c1 + sqrt_disc) / 2
            else:
                r1 = (-c1 - sqrt_disc) / 2
            r2 = c2 / r1
            return [r1, r2]
        else:
            return []

    def roots_c(self):
        # second order sections (2x roots either real or complex conj)

        # complex form is always numerically stable
        c1, c2 = self.effective_params()
        # a = c2, b = c1, c = 1
        disc = c1 * c1 - 4 * c2
        if disc >= 0:
            sqrt_disc = disc ** 0.5
            return []
        else:
            sqrt_disc = (-disc) ** 0.5
            r1 = (-c1 + sqrt_disc * 1j) / 2
            return [r1]


class CodingSOSMirror(CodingSOS):
    unstable = False
