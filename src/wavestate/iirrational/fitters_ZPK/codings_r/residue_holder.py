#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import numpy as np

from .base import (
    CodingType,
    Ipi,
    #I2pi
)
from ... import TFmath
#import scipy.linalg


class CodingResidues(CodingType):
    N_parameters = None

    def setup(self, coding_K, codings_r):
        self.coding_K = coding_K
        self.codings_r = codings_r
        self.codings_full = [coding_K, codings_r]

        N_parameters = 0
        for coding in self.codings_full:
            N_parameters += coding.N_parameters

        self.N_parameters = N_parameters
        return

    #def update_ZPK(self, ZPK):
    #    res_list, ZPKsum = TFmath.ZPK2residues()
    #    return

    def update(self, *params):
        idx_p = 0
        for coding in self.codings_full:
            subP = params[idx_p:idx_p + coding.N_parameters]
            coding.update(*subP)
            idx_p += coding.N_parameters
        return

    def reduce(self):
        params = []
        for coding in self.codings_full:
            params.extend(
                coding.reduce()
            )
        return params

    def option_set(self, **kwargs):
        for coding in self.codings_full:
            coding.option_set(**kwargs)
        return

    def transfer(self, sys):
        val = 0
        for coding in self.codings_full:
            subT = coding.transfer(sys)
            val = val + subT
        return subT

    def derivative(self, sys):
        derivs = []
        for coding in self.codings_full:
            subD = coding.derivative(sys)
            derivs.extend(subD)
        return derivs

    def ZPK(self, sys):
        residue_list = []
        for c_r in self.codings_r:
            residue_list.extend(
                c_r.residues(sys)
            )

        ZPKsum = self.coding_K.ZPK(sys)

        return TFmath.residues2ZPK(residue_list, ZPKsum)


