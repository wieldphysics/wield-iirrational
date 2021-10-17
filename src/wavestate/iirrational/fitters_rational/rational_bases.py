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
from declarative import (
    DepBunch,
    #Bunch,
    depB_property,
    NOARG,
)

from .. import representations

class DataFilterBase(DepBunch):
    RBalgo = representations.RBAlgorithms()
    root_constraint = RBalgo.root_constraints.mirror_real

    def __build__(
        self,
        _args  = None,
        parent = NOARG,
        ZPKrep = NOARG,
        **kwargs
    ):
        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")

        if parent is not NOARG and ZPKrep is NOARG:
            ZPKrep = parent.ZPKrep

        if ZPKrep is not NOARG:
            kwargs['parent'] = ZPKrep
        ZPKrep = representations.ZPKwData(**kwargs)

        super(DataFilterBase, self).__build__()
        self.ZPKrep = ZPKrep

    @depB_property
    def F_max_Hz(self):
        return np.max(self.F_Hz)

    @depB_property
    def X_grid(self):
        return self.phi_rep(self.F_Hz)

    def xfer_eval(self, F_Hz):
        #TODO, use ZPKrep for this...
        alt_filter = self.__class__(
            F_Hz   = F_Hz,
            W      = None,
            data   = None,
            parent = self,
        )
        #providing the numerical scale will aid its precision
        #note! it should be computed in a more numerically stable way!
        #TODO use a ZPK computation for this
        return alt_filter.xfer_fit

    @depB_property
    def ZPKrep(self, rep = NOARG):
        if rep is NOARG:
            Z, lnGZ     = self.phi_rep_native_lnG(self.zeros)
            P, lnGP     = self.phi_rep_native_lnG(self.poles)
            Zov, lnGZov = self.phi_rep_native_lnG(self.zeros_overlay)
            Pov, lnGPov = self.phi_rep_native_lnG(self.poles_overlay)
            return representations.ZPKwData(
                data          = self.data,
                W             = self.W,
                F_Hz          = self.F_Hz,
                zeros         = Z,
                poles         = P,
                zeros_overlay = Zov,
                poles_overlay = Pov,
                gain          = self.gain * np.exp(-(lnGZ + lnGZov) + (lnGP + lnGPov)),
                delay_s       = self.delay_s,
                F_nyquist_Hz  = self.F_nyquist_Hz,
            )
        else:
            self.data         = rep.data
            self.W            = rep.W
            self.F_Hz         = rep.F_Hz
            self.delay_s      = rep.delay_s
            self.F_nyquist_Hz = rep.F_nyquist_Hz

            #must be assigned
            Z, lnGZ     = self.phi_native_rep_lnG(rep.zeros)
            P, lnGP     = self.phi_native_rep_lnG(rep.poles)
            Zov, lnGZov = self.phi_native_rep_lnG(rep.zeros_overlay)
            Pov, lnGPov = self.phi_native_rep_lnG(rep.poles_overlay)
            self.zeros = Z
            self.poles = P
            self.zeros_overlay = Zov
            self.poles_overlay = Pov
            self.gain = rep.gain * np.exp(-(lnGZ + lnGZov) + (lnGP + lnGPov))

            self.dependencies('gain')
            self.dependencies('zeros')
            self.dependencies('poles')
            self.dependencies('zeros_overlay')
            self.dependencies('poles_overlay')
            self.dependencies('data')
            self.dependencies('W')
            self.dependencies('delay_s')
            self.dependencies('F_Hz')
            self.dependencies('F_nyquist_Hz')
            return rep

    def phi_coeff_norm(self, pvec):
        return pvec[-1] * self.phi_gain_adj_pvec(pvec)

    def phi_gain_adj_pvec(self, pvec):
        return 1

    def phi_gain_adj_roots(self, rB):
        return 1

    def phi_rep(self, F_Hz):
        """
        function mapping data frequencies from frequency domain to fit domain
        """
        return 1j * F_Hz

    def phi_rep_native_lnG(self, rB):
        """
        function mapping roots to native domain (S) from fit representation domain
        outputs both the remapped roots, as well as the lnG adjustment
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB, 0

    def phi_native_rep_lnG(self, rB):
        """
        function mapping roots from native domain (S) to fit representation domain
        outputs both the remapped roots, as well as the lnG adjustment
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB, 0

    def phi_rep_Snative(self, rB):
        """
        Maps the internal representation to Snative, the output RootBunch can
        have a specialty constraint added.
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB

    def phi_Snative_rep(self, rB):
        """
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB

    def _phasing_roots(self, roots):
        return 1

    def _phasing_vec(self, vec):
        return 1

    def _phasing_afit(self):
        return 1

    def _phasing_bfit(self):
        return 1

    @depB_property
    def order(self):
        return max(
            len(self.poles) + len(self.poles_overlay),
            len(self.zeros) + len(self.zeros_overlay),
        )

    @depB_property
    def order_sos(self):
        return (max(
            len(self.poles) + len(self.poles_overlay),
            len(self.zeros) + len(self.zeros_overlay),
        ) + 1)//2

    @depB_property
    def order_relative(self):
        return (
            (len(self.zeros) + len(self.zeros_overlay))
            - (len(self.poles) + len(self.poles_overlay))
        )



class ZFilterBase(DataFilterBase):
    @depB_property
    def F_nyquist_Hz(self, val = NOARG):
        if val is NOARG:
            raise RuntimeError("Must specify F_nyquist_Hz before using class")
        if val is None:
            raise RuntimeError("F_nyquist_Hz must be positive for Z filters")
        return val

    @depB_property
    def Xn_grid(self):
        return self.phi_rep(-self.F_Hz)

    @depB_property
    def ZPKz(self):
        rep = self.ZPKrep
        return representations.ZPKCalc(
            tuple(rep.zeros.fullplane) + tuple(rep.zeros_overlay.fullplane),
            tuple(rep.poles.fullplane) + tuple(rep.poles_overlay.fullplane),
            self.gain,
            F_nyquist_Hz = self.F_nyquist_Hz,
        )

    def phi_gain_adj_roots(self, roots):
        return 1

    def phi_gain_adj_pvec(self, pvec):
        return 1

    def phi_rep(self, F_Hz):
        """
        function mapping roots from frequency domain to fit domain
        """
        return np.exp(1j * np.pi * F_Hz / self.F_nyquist_Hz)

    def phi_rep_native_lnG(self, rB):
        """
        function mapping roots to native domain (S) from fit representation domain
        outputs both the remapped roots, as well as the lnG adjustment
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB, 0

    def phi_native_rep_lnG(self, rB):
        """
        function mapping roots from native domain (S) to fit representation domain
        outputs both the remapped roots, as well as the lnG adjustment
        """
        rB = self.RBalgo.expect(rB, self.RBalgo.root_constraints.mirror_real)
        return rB, 0

    def phi_rep_Snative(self, rB):
        """
        Maps the internal representation to Snative, the output RootBunch can
        have a specialty constraint added.
        """
        rB.c
        return rB

    def phi_Snative_rep(self, rB):
        """
        """
        return rB


class SFilterBase(DataFilterBase):
    def __init__(
        self,
        _args  = None,
        F_nyquist_Hz = None,
        **kwargs
    ):
        super(SFilterBase, self).__init__(
            _args,
            F_nyquist_Hz = F_nyquist_Hz,
            **kwargs
        )

    @depB_property
    def F_nyquist_Hz(self, val = NOARG):
        if val is NOARG:
            return None
        if val is not None:
            raise RuntimeError("F_nyquist_Hz must None for S filters")
        return val

    @depB_property
    def ZPKsf(self):
        rep = self.ZPKrep
        return representations.ZPKCalc(
            tuple(rep.zeros.fullplane) + tuple(rep.zeros_overlay.fullplane),
            tuple(rep.poles.fullplane) + tuple(rep.poles_overlay.fullplane),
            rep.gain,
            F_nyquist_Hz = rep.F_nyquist_Hz,
        )
