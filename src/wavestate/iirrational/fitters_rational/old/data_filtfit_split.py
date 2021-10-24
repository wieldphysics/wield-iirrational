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
from wavestate.bunch.depbunch import DepBunch, Bunch

from ..roots_bin import roots_bin_type
from ..polynomials import poly_group


def abs_sq(x):
    return x.real ** 2 + x.imag ** 2


class DataFiltFitSplitBase(DepBunch):
    def __build__(
        self,
        _args=None,
        F_Hz=None,
        data=None,
        W=None,
        Fz_nyquist_Hz=None,
        Fp_nyquist_Hz=None,
        F_ref_Hz=None,
        delay_s=None,
        gain=None,
        **kwargs
    ):
        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")
        if F_Hz is None:
            raise RuntimeError("Must Specify F_Hz")
        if data is None:
            raise RuntimeError("Must Specify data")
        if W is None:
            W = 1
        if Fp_nyquist_Hz is None:
            raise RuntimeError("Must Specify Fp_nyquist_Hz")
        if Fz_nyquist_Hz is None:
            raise RuntimeError("Must Specify Fz_nyquist_Hz")
        if F_ref_Hz is None:
            F_ref_Hz = 0
        if delay_s is None:
            delay_s = 0
        if gain is None:
            gain = 1

        if data is not None and len(np.asarray(data).flatten()) == 0:
            data = None
        if W is not None and len(np.asarray(W).flatten()) == 0:
            W = None

        if data is not None:
            if W is None:
                W = 1
            F_Hz, data, W = np.broadcast_arrays(F_Hz, data, W)
        self.F_ref_Hz = F_ref_Hz
        self.Fz_nyquist_Hz = Fz_nyquist_Hz
        self.Fp_nyquist_Hz = Fp_nyquist_Hz
        self.F_Hz = F_Hz
        self.data = data
        self.W = W
        self.delay_s = delay_s
        self.gain = gain

        super(DataFiltFitSplitBase, self).__build__(**kwargs)

        @self.deco_generator
        def poly(self):
            return poly_group(self.Fp_nyquist_Hz)

        @self.deco_generator
        def poles_split(self):
            p_r, p_c, p_u, pol = roots_bin_type(
                self.poles, policy="strict", F_nyquist_Hz=self.Fp_nyquist_Hz
            )
            if len(p_u) > 0:
                raise RuntimeError(
                    "Got unorganized poles (missing conjugates). Not Prepared!"
                )
            return p_r, p_c

        @self.deco_generator
        def zeros_split(self):
            z_r, z_c, z_u, pol = roots_bin_type(
                self.zeros, policy="strict", F_nyquist_Hz=self.Fz_nyquist_Hz
            )
            if len(z_u) > 0:
                raise RuntimeError(
                    "Got unorganized zeros (missing conjugates). Not Prepared!"
                )
            return z_r, z_c

        @self.deco_generator
        def residuals(self):
            if self.data is None:
                raise RuntimeError("Can't generate residuals as data was not Specified")
            debias_reweight = 1 / (0.001 + self.W ** 2)
            retB = Bunch()
            retB.resP = self.W * (self.xfer_fit / self.data - 1)
            retB.resZ = self.W * (self.data / self.xfer_fit - 1)
            retB.total = np.sum(
                (abs(retB.resP) ** 2 + abs(retB.resZ * debias_reweight) ** 2)
                / (1 + debias_reweight)
            ) / (2 * len(self.data))
            return retB

        @self.deco_generator
        def F_max_Hz(self):
            return np.max(self.F_Hz)

        @self.deco_generator
        def Xp_grid_P(self):
            return np.exp(1j * np.pi * self.F_Hz / self.Fp_nyquist_Hz)

        @self.deco_generator
        def Xp_ref_P(self):
            return np.exp(1j * np.pi * self.F_ref_Hz / self.Fp_nyquist_Hz)

        @self.deco_generator
        def Xn_grid_P(self):
            return np.exp(-1j * np.pi * self.F_Hz / self.Fp_nyquist_Hz)

        @self.deco_generator
        def Xn_ref_P(self):
            return np.exp(-1j * np.pi * self.F_ref_Hz / self.Fp_nyquist_Hz)

        @self.deco_generator
        def Xp_grid_Z(self):
            return np.exp(1j * np.pi * self.F_Hz / self.Fz_nyquist_Hz)

        @self.deco_generator
        def Xp_ref_Z(self):
            return np.exp(1j * np.pi * self.F_ref_Hz / self.Fz_nyquist_Hz)

        @self.deco_generator
        def Xn_grid_Z(self):
            return np.exp(-1j * np.pi * self.F_Hz / self.Fz_nyquist_Hz)

        @self.deco_generator
        def Xn_ref_Z(self):
            return np.exp(-1j * np.pi * self.F_ref_Hz / self.Fz_nyquist_Hz)

    @property
    def ZPK(self):
        # cheap to compute, no need to autodelete
        return (self.zeros, self.poles, self.gain)
