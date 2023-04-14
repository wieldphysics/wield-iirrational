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
from wavestate.bunch.depbunch import (
    DepBunch,
    depB_property,
    NOARG,
)

from .root_bunch import RBAlgorithms
from .zpktf import ZPKTF


class ZPKwData(DepBunch):
    RBalgo = RBAlgorithms()
    root_constraint = RBalgo.root_constraints.mirror_real
    phase_missing = False
    residuals_average = 0

    def __build__(
        self,
        _args=None,
        data=NOARG,
        F_Hz=NOARG,
        W=NOARG,
        delay_s=NOARG,
        zeros=NOARG,
        poles=NOARG,
        gain=NOARG,
        zeros_overlay=NOARG,
        poles_overlay=NOARG,
        ZPK=NOARG,
        ZPK_overlay=NOARG,
        F_nyquist_Hz=NOARG,
        parent=NOARG,
        residuals_log_im_scale=NOARG,
        ZPKrep=NOARG,
    ):

        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")

        if zeros is NOARG and ZPK is not NOARG:
            zeros = ZPK[0]

        if poles is NOARG and ZPK is not NOARG:
            poles = ZPK[1]

        if zeros_overlay is NOARG and ZPK_overlay is not NOARG:
            zeros_overlay = ZPK_overlay[0]

        if poles_overlay is NOARG and ZPK_overlay is not NOARG:
            poles_overlay = ZPK_overlay[1]

        if gain is NOARG and ZPK is not NOARG:
            gain = ZPK[2]

        if ZPK_overlay is not NOARG:
            gain2 = ZPK_overlay[2]
        else:
            gain2 = 1

        if ZPKrep is not NOARG:
            parent = ZPKrep

        if parent is not NOARG:
            if data is NOARG:
                data = parent.data
            if F_Hz is NOARG:
                F_Hz = parent.F_Hz
            if W is NOARG:
                W = parent.W
            if delay_s is NOARG:
                delay_s = parent.delay_s
            if zeros is NOARG:
                zeros = parent.zeros
            if poles is NOARG:
                poles = parent.poles
            if gain is NOARG:
                gain = parent.gain
            if zeros_overlay is NOARG:
                zeros_overlay = parent.zeros_overlay
            if poles_overlay is NOARG:
                poles_overlay = parent.poles_overlay
            if F_nyquist_Hz is NOARG:
                F_nyquist_Hz = parent.F_nyquist_Hz
            if residuals_log_im_scale is NOARG:
                residuals_log_im_scale = parent.residuals_log_im_scale

        if zeros is NOARG:
            zeros = ()

        if poles is NOARG:
            poles = ()

        if zeros_overlay is NOARG:
            zeros_overlay = ()

        if poles_overlay is NOARG:
            poles_overlay = ()

        if gain is NOARG:
            gain = 1

        if F_Hz is NOARG:
            raise RuntimeError("Must Specify F_Hz")
        if data is NOARG:
            raise RuntimeError("Must Specify data")
        if F_nyquist_Hz is NOARG:
            raise RuntimeError("Must Specify F_nyquist_Hz")
        if W is NOARG:
            if data is None:
                W = None
            else:
                W = 1
        if residuals_log_im_scale is NOARG:
            if data is None:
                residuals_log_im_scale = None
            else:
                residuals_log_im_scale = 1

        if delay_s is NOARG:
            delay_s = 0

        if data is not None:
            if W is None:
                W = 1
            F_Hz, data, W, residuals_log_im_scale = np.broadcast_arrays(
                F_Hz, data, W, residuals_log_im_scale
            )

        self.F_Hz = F_Hz
        self.data = data
        self.W = W
        self.delay_s = delay_s

        self.zeros = zeros
        self.poles = poles
        self.zeros_overlay = zeros_overlay
        self.poles_overlay = poles_overlay
        self.gain = gain * gain2

        self.F_nyquist_Hz = F_nyquist_Hz

        self.residuals_log_im_scale = residuals_log_im_scale

    @depB_property
    def order(self):
        return max(
            len(self.poles) + len(self.poles_overlay),
            len(self.zeros) + len(self.zeros_overlay),
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
    def F_nyquist_Hz(self, val):
        if val is None:
            return val
        assert val > 0
        return val

    @depB_property
    def zeros(self, val):
        val = self.RBalgo.expect_atleast(val, constraint=self.root_constraint)
        return val

    @depB_property
    def poles(self, val):
        val = self.RBalgo.expect_atleast(val, constraint=self.root_constraint)
        return val

    @depB_property
    def zeros_overlay(self, val):
        val = self.RBalgo.expect_atleast(val, constraint=self.root_constraint)
        return val

    @depB_property
    def poles_overlay(self, val):
        val = self.RBalgo.expect_atleast(val, constraint=self.root_constraint)
        return val

    @depB_property
    def X_grid(self):
        if self.F_nyquist_Hz is None:
            return 1j * self.F_Hz
        else:
            # use Z^-1
            return np.exp(1j * np.pi * self.F_Hz / self.F_nyquist_Hz)

    @depB_property
    def ZPKrep(self):
        return self

    @depB_property
    def ZPK(self):
        return ZPKTF(
            zeros=self.zeros * self.zeros_overlay,
            poles=self.poles * self.poles_overlay,
            gain=self.gain,
            F_nyquist_Hz=self.F_nyquist_Hz,
        )

    @depB_property
    def ZPKfit(self):
        return ZPKTF(
            (
                self.zeros,
                self.poles,
                self.gain,
            ),
            F_nyquist_Hz=self.F_nyquist_Hz,
        )

    @depB_property
    def ZPKoverlay(self):
        return ZPKTF(
            (
                self.zeros_overlay,
                self.poles_overlay,
                1,
            ),
            F_nyquist_Hz=self.F_nyquist_Hz,
        )

    @depB_property
    def xfer_fit(self):
        # TODO must add pole-zero rephasing for Z filters

        # TODO, make this name consistent in all classes
        h, lnG = self.poles_overlay.val_lnG(self.X_grid)
        h, lnG = self.poles.val_lnG(self.X_grid, h=h, lnG=lnG)

        h, lnG = self.zeros_overlay.val_lnG(self.X_grid, h=1 / h, lnG=-lnG)
        h, lnG = self.zeros.val_lnG(self.X_grid, h=h, lnG=lnG)
        if self.delay_s != 0:
            h = h * np.exp(-2j * np.pi * self.delay_s * self.F_Hz)
        return h * (np.exp(lnG) * self.gain)

    def xfer_eval(self, F_Hz):
        # TODO must add pole-zero rephasing for Z filters

        # TODO, make this name consistent in all classes
        if self.F_nyquist_Hz is None:
            X_grid = 1j * F_Hz
        else:
            # use Z^-1
            X_grid = np.exp(1j * np.pi * F_Hz / self.F_nyquist_Hz)
        h, lnG = self.poles_overlay.val_lnG(X_grid)
        h, lnG = self.poles.val_lnG(X_grid, h=h, lnG=lnG)

        h, lnG = self.zeros_overlay.val_lnG(X_grid, h=1 / h, lnG=-lnG)
        h, lnG = self.zeros.val_lnG(X_grid, h=h, lnG=lnG)
        if self.delay_s != 0:
            h = h * np.exp(-2j * np.pi * self.delay_s * F_Hz)
        return h * (np.exp(lnG) * self.gain)

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
