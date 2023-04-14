#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import numpy as np
from wield.bunch.depbunch import depB_property, NOARG

from .. import representations
from ..representations.polynomials import standard

from .rational_algorithms import PolyFilterAlgorithms
from .rational_attributes import PolyFit
from .rational_bases import ZFilterBase


class RootConstraints(representations.RootConstraints):
    branch_separated = frozenset("branch_separated")


root_constraints = RootConstraints()


# must be in this order for the ZPKz to be set in the correct order
class RationalDiscFilter(ZFilterBase, PolyFit, PolyFilterAlgorithms):
    phase_missing = False
    poly = standard
    stablize_pos_method = "ignore"
    stablize_neg_method = "ignore"

    @depB_property
    def zeros_phase_adjust(self, val=NOARG):
        if val is NOARG:
            return None
        else:
            return val

    @depB_property
    def poles_phase_adjust(self, val=NOARG):
        if val is NOARG:
            return None
        else:
            return val

    @depB_property
    def zeros_phasing(self):
        if self.zeros_phase_adjust is None:
            return self.Xn_grid ** self.nzeros
        elif self.zeros_phase_adjust < 0:
            return self.Xn_grid ** (self.nzeros + self.zeros_phase_adjust)
        else:
            return self.Xn_grid ** (self.zeros_phase_adjust)
        return

    @depB_property
    def poles_phasing(self):
        if self.poles_phase_adjust is None:
            return self.Xn_grid ** self.npoles
        elif self.poles_phase_adjust < 0:
            return self.Xn_grid ** (self.npoles + self.poles_phase_adjust)
        else:
            return self.Xn_grid ** (self.poles_phase_adjust)
        return

    def _phasing_roots(self, roots):
        return self.Xn_grid ** len(roots)

    def _phasing_vec(self, vec):
        return self.Xn_grid ** (len(vec) - 1)

    def _phasing_afit(self):
        return self.poles_phasing

    def _phasing_bfit(self):
        return self.zeros_phasing

    def root_stabilize(
        self,
        rB,
        method=None,
        real_neg=None,
        real_pos=None,
    ):
        if method is None:
            return rB

        if real_neg is None:
            real_neg = self.stablize_neg_method
        if real_pos is None:
            real_pos = self.stablize_pos_method

        r_r = np.copy(rB.r)
        r_c = np.copy(rB.c)

        if real_neg == "use":
            select_r = r_r < -1
        elif real_neg == "ignore":
            select_r = np.zeros(r_r.shape, dtype=bool)
        else:
            raise RuntimeError("Bad Argument")

        if real_pos == "use":
            select_r = r_r > 1 | select_r
        elif real_pos == "ignore":
            pass
        else:
            raise RuntimeError("Bad Argument")

        select_c = abs(r_c) > 1
        if method == "flip":
            r_r[select_r] = 1 / r_r[select_r]
            r_c[select_c] = 1 / r_c[select_c].conjugate()
        elif method == "remove":
            r_r = r_r[~select_r]
            r_c = r_c[~select_c]
        else:
            raise RuntimeError("Bad Argument")

        rB = representations.RootBunch(
            constraint=self.root_constraint,
            r=r_r,
            c=r_c,
        )
        return rB

    def clear_unstable_poles(self):
        new_poles = []
        for r in self.poles:
            if abs(r) < 1:
                new_poles.append(r)
        self.poles = new_poles

    def clear_bigdelay_zeros(self):
        new_zeros = []
        for r in self.zeros:
            if abs(r) < 1:
                new_zeros.append(r)
        self.zeros = new_zeros

    def phi_rep_Snative(self, rB):
        """
        Maps the internal representation to Snative, the output RootBunch can
        have a specialty constraint added.
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        # TODO
        # ok to modify?
        rB.constraint = rB.constraint | root_constraints.branch_separated
        select_pos = rB.r > 0
        rB.nyquist_branch = (-rB.r[~select_pos] - 1) * self.F_nyquist_Hz
        rB.r = (rB.r[select_pos] - 1) * self.F_nyquist_Hz
        real = (abs(rB.c) - 1) * self.F_nyquist_Hz
        imag = np.angle(rB.c) / np.pi * self.F_nyquist_Hz
        rB.c = real + 1j * imag
        return rB

    def phi_Snative_rep(self, rB):
        """ """
        assert root_constraints.mirror_real == (
            rB.constraint - root_constraints.branch_separated
        )
        rB.c = (1 + rB.c.real / self.F_nyquist_Hz) * np.exp(
            1j * rB.c.imag * np.pi / self.F_nyquist_Hz
        )

        if root_constraints.branch_separated <= rB.constraint:
            rB.r = np.concatenate(
                [
                    (rB.r / self.F_nyquist_Hz + 1),
                    -(rB.nyquist_branch / self.F_nyquist_Hz + 1),
                ]
            )
        else:
            rB.r = rB.r / self.F_nyquist_Hz + 1

        rB.constraint = rB.constraint - root_constraints.branch_separated
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB
