# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np

from ..representations.polynomials import (
    chebychev,
    poly_constraints
)

#from .. import TFmath
from .rational_attributes import PolyFit
from .rational_bases import SFilterBase
from .rational_algorithms import PolyFilterAlgorithmsIm


#must be in this order for the ZPKsf to be set in the proper order
class ChebychevFilter(
    SFilterBase,
    PolyFit,
    PolyFilterAlgorithmsIm
):
    phase_missing = False
    root_constraint = SFilterBase.RBalgo.root_constraints.mirror_imag
    poly = chebychev

    def phi_gain_adj_roots(self, roots):
        if len(roots) % 2:
            return 1j
        else:
            return 1

    def phi_gain_adj_pvec(self, pvec):
        if len(pvec) % 2:
            return 1
        else:
            return 1j

    def phi_rep(self, F_Hz):
        """
        function mapping roots from frequency domain to fit domain
        """
        return F_Hz / self.F_max_Hz

    def phi_rep_native_lnG(self, rB):
        """
        function mapping roots to native domain (S) from fit representation domain
        outputs both the remapped roots, as well as the lnG adjustment
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_imag)
        lnG = len(rB) * np.log(self.F_max_Hz)
        return rB * (1j * self.F_max_Hz), lnG

    def phi_native_rep_lnG(self, rB):
        """
        function mapping roots from native domain (S) to fit representation domain
        outputs both the remapped roots, as well as the lnG adjustment
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        lnG = len(rB) * np.log(self.F_max_Hz)
        return rB * (-1j / self.F_max_Hz), -lnG

    def phi_rep_Snative(self, rB):
        """
        Maps the internal representation to Snative, the output RootBunch can
        have a specialty constraint added.
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_imag)
        return rB * (1j * self.F_max_Hz)

    def phi_Snative_rep(self, rB):
        """
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB * (-1j / self.F_max_Hz)

