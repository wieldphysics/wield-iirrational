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
from wield.bunch.depbunch import DepBunch, depB_property
from wield import declarative
from wield.bunch import Bunch


from .. import representations


class DataFiltFitBase(DepBunch):
    def __build__(
        self,
        _args=None,
        F_Hz=None,
        data=None,
        W=None,
        parent=None,
        codings=None,
        F_cutoff_Hz=None,
        **kwargs
    ):
        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")
        if parent is not None:
            if F_Hz is None:
                F_Hz = parent.F_Hz
            if data is None:
                data = parent.data
            if W is None:
                W = parent.W
            if codings is None:
                codings = parent.codings
            if F_cutoff_Hz is None:
                F_cutoff_Hz = parent.F_cutoff_Hz
        else:
            if F_Hz is None:
                raise RuntimeError("Must Specify F_Hz")
            if data is None:
                raise RuntimeError("Must Specify data")
            if W is None:
                W = 1

        if data is not None and len(np.asarray(data).flatten()) == 0:
            data = None
        if W is not None and len(np.asarray(W).flatten()) == 0:
            W = None

        if data is not None:
            if W is None:
                W = 1
            F_Hz, data, W = np.broadcast_arrays(F_Hz, data, W)
        self.F_Hz = F_Hz
        self.data = data
        self.W = W
        self.codings = codings

        if F_cutoff_Hz is not None:
            self.F_cutoff_Hz = F_cutoff_Hz

        super(DataFiltFitBase, self).__build__(**kwargs)

    @depB_property
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

    @depB_property
    def F_max_Hz(self):
        return np.max(self.F_Hz)

    @depB_property
    def F_cutoff_Hz(self):
        return 1.00 * np.max(self.F_Hz)

    @depB_property
    def order(self):
        return max(
            len(self.poles) + len(self.poles_overlay),
            len(self.zeros) + len(self.zeros_overlay),
        )

    @depB_property
    def order_total(self):
        return (
            len(self.poles)
            + len(self.poles_overlay)
            + len(self.zeros)
            + len(self.zeros_overlay)
        )

    @depB_property
    def order_sos(self):
        return (
            max(
                len(self.poles) + len(self.poles_overlay),
                len(self.zeros) + len(self.zeros_overlay),
            )
            + 1
        ) // 2

    @depB_property
    def order_relative(self):
        return (len(self.zeros) + len(self.zeros_overlay)) - (
            len(self.poles) + len(self.poles_overlay)
        )

    @depB_property
    def ZPKrep(self):
        # cheap to compute, no need to autodelete
        self.codings_revision
        return representations.ZPKwData(
            data=self.data,
            W=self.W,
            F_Hz=self.F_Hz,
            zeros=self.zeros,
            poles=self.poles,
            zeros_overlay=self.zeros_overlay,
            poles_overlay=self.poles_overlay,
            gain=self.gain,
            delay_s=self.delay_s,
            F_nyquist_Hz=self.F_nyquist_Hz,
        )


class DataFiltFitZ(DataFiltFitBase):
    def __build__(
        self, _args=None, F_nyquist_Hz=declarative.NOARG, parent=None, **kwargs
    ):
        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")
        if parent is not None:
            if F_nyquist_Hz is declarative.NOARG:
                F_nyquist_Hz = parent.F_nyquist_Hz
        else:
            if F_nyquist_Hz is declarative.NOARG:
                raise RuntimeError("Must Specify F_nyquist_Hz")

        if F_nyquist_Hz is None:
            raise RuntimeError("This type cannot support F_nyquist_Hz = None")

        self.F_nyquist_Hz = F_nyquist_Hz

        super(DataFiltFitZ, self).__build__(parent=parent, **kwargs)

    @depB_property
    def Xzp_grid(self):
        return np.exp(1j * np.pi * self.F_Hz / self.F_nyquist_Hz)

    @depB_property
    def Xzp_grid_sq(self):
        return self.Xzp_grid * self.Xzp_grid

    @depB_property
    def Xzn_grid(self):
        return np.exp(-1j * np.pi * self.F_Hz / self.F_nyquist_Hz)

    @depB_property
    def Xzn_grid_sq(self):
        return self.Xzn_grid * self.Xzn_grid

    @depB_property
    def Xex_grid(self):
        return self.Xzp_grid


class DataFiltFitSf(DataFiltFitBase):
    def __build__(
        self, _args=None, F_nyquist_Hz=declarative.NOARG, parent=None, **kwargs
    ):
        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")
        if parent is not None:
            if F_nyquist_Hz is declarative.NOARG:
                F_nyquist_Hz = parent.F_nyquist_Hz
        else:
            if F_nyquist_Hz is declarative.NOARG:
                F_nyquist_Hz = None

        if F_nyquist_Hz is not None:
            raise RuntimeError("This type cannot support F_nyquist_Hz != None")

        self.F_nyquist_Hz = F_nyquist_Hz

        super(DataFiltFitSf, self).__build__(parent=parent, **kwargs)

    @depB_property
    def Xsf_grid(self):
        return 1j * self.F_Hz

    @depB_property
    def Xsf_grid_sq(self):
        return self.Xsf_grid * self.Xsf_grid

    @depB_property
    def Xex_grid(self):
        return self.Xsf_grid
