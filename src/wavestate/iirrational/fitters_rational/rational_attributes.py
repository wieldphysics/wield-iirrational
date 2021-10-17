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
import warnings
import declarative
from declarative import (
    depB_property,
    NOARG
)
import collections

from .. import TFmath

from .rational_bases import DataFilterBase


lnGtup = collections.namedtuple('lnGtup', ['sign', 'lnG'])
lnGxfer = collections.namedtuple('lnGxfer', ['xfer', 'lnG'])
lnGVand = collections.namedtuple('lnGVand', ['vander', 'lnG'])


class PolyFit(DataFilterBase):
    phase_missing = False
    #must fill
    poly = None
    roots_split_strict = False
    stabilize_method = 'flip'

    def __build__(
        self,
        _args     = None,
        npoles    = NOARG,
        nzeros    = NOARG,
        mindelay  = False,
        stabilize = False,
        parent    = NOARG,
        **kwargs
    ):
        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")
        if parent is not NOARG:
            if npoles is NOARG:
                npoles = parent.npoles
            if nzeros is NOARG:
                nzeros = parent.nzeros
            if mindelay is NOARG:
                npoles = parent.mindelay
            if nzeros is NOARG:
                stabilize = parent.stabilize

        self.npoles    = npoles
        self.nzeros    = nzeros
        self.mindelay  = mindelay
        self.stabilize = stabilize

        super(PolyFit, self).__build__(
            parent = parent,
            **kwargs
        )
        #use these for now
        #self.lnG = 0, 1
        #self.gain = 1

        ##first to prevent lookups that shouldn't happen
        #self.zeros     = []
        #self.poles     = []
        #print(self.zeros)
        #print(self.poles)
        #self.a_vec
        #self.b_vec

    @depB_property
    def lnG(self, val = NOARG):
        if val is NOARG:
            raise RuntimeError("shouldn't compute")
            return self.b_vec[-1]/self.a_vec[-1]
        else:
            sgn, lnG = val
            if not np.isfinite(lnG):
                raise RuntimeError("Gain must be finite")
            return lnGtup(sgn, lnG)

    @depB_property
    def gain(self, val = NOARG):
        if val is NOARG:
            return self.lnG.sign * np.exp(self.lnG.lnG)
        else:
            if val == 0:
                raise RuntimeError("Gain 0 not allowed")
            if val.imag != 0:
                raise RuntimeError("Only Real Gain allowed")
            elif not np.isfinite(val):
                raise RuntimeError("Gain Must be finite")
            elif val > 0:
                self.lnG = 1, np.log(val)
            else:
                self.lnG = -1, np.log(-val)
            self.dependencies('lnG')
            return val

    @depB_property
    def a_vec(self, pvec = NOARG):
        if pvec is NOARG:
            vec = self.poly_fromroots(self.poles)
            norm = self.phi_coeff_norm(vec)
            if norm != 1:
                warnings.warn("poly normalization unexpected, {0}".format(norm))
            return vec
        else:
            pvec = np.asarray(pvec)
            #adjust gain numbers and canonicalize
            signG, lnG = self.get_raw('lnG')
            cG = self.phi_coeff_norm(pvec)
            self.lnG = signG * np.sign(cG), lnG - np.log(abs(cG))
            pvec = pvec / cG

            if self.get_raw('stabilize'):
                poles = self.poly.roots_rB(pvec, self.poly_constraints)
                #only prevent using this vector if needed, stability is helped if we don't go convert between roots and vectors too much
                #actually set it so that we get the stabilize action
                self.poles = poles
                #don't store any value as it should be generated from the poles if needed
                return self.TAG_NO_SET

            self.dependencies('poles')
            return pvec

    @depB_property
    def b_vec(self, pvec = NOARG):
        if pvec is NOARG:
            vec = self.poly_fromroots(self.zeros)
            norm = self.phi_coeff_norm(vec)
            if norm != 1:
                warnings.warn("poly normalization unexpected, {0}".format(norm))
            return vec
        else:
            pvec = np.asarray(pvec)

            signG, lnG = self.get_raw('lnG')
            cG = self.phi_coeff_norm(pvec)
            self.lnG = signG * np.sign(cG), lnG + np.log(abs(cG))
            pvec = pvec / cG

            if self.get_raw('mindelay'):
                zeros = self.poly.roots_rB(pvec, self.poly_constraints)
                #only prevent using this vector if needed, stability is helped if we don't go convert between roots and vectors too much
                #actually set it so that we get the stabilize action
                self.zeros = zeros
                return self.TAG_NO_SET
            #since they had to be computed anyway for the gain adjustment
            self.dependencies('zeros')
            return pvec

    @depB_property
    def zeros(self, val = NOARG):
        if val is NOARG:
            rB = self.poly.roots_rB(self.b_vec, self.poly_constraints)
            return self.RBalgo.expect(rB, self.root_constraint)
        else:
            #get previous
            rB = self.RBalgo.expect(val, self.root_constraint)
            if self.get_raw('mindelay'):
                rB = self.root_stabilize(
                    rB,
                    method = self.stabilize_method,
                )
            self.dependencies('b_vec')
            return rB

    @depB_property
    def poles(self, val = NOARG):
        if val is NOARG:
            rB = self.poly.roots_rB(self.a_vec, self.poly_constraints)
            return self.RBalgo.expect(rB, self.root_constraint)
        else:
            rB = self.RBalgo.expect(val, self.root_constraint)
            if self.get_raw('stabilize'):
                rB = self.root_stabilize(
                    rB,
                    method = self.stabilize_method,
                )
            self.dependencies('a_vec')
            return rB

    @depB_property
    def zeros_overlay(self, val = NOARG):
        if val is NOARG:
            return self.RBalgo.expect([], self.root_constraint)
        else:
            rB = self.RBalgo.expect(val, self.root_constraint)
            return rB

    @depB_property
    def poles_overlay(self, val = NOARG):
        if val is NOARG:
            return self.RBalgo.expect([], self.root_constraint)
        else:
            rB = self.RBalgo.expect(val, self.root_constraint)
            return rB

    @depB_property
    def V_b(self):
        return lnGVand(*self.poly.vander_lnG(self.X_grid, self.nzeros))

    @depB_property
    def V_a(self):
        return lnGVand(*self.poly.vander_lnG(self.X_grid, self.npoles))

    @depB_property
    def h_delay(self):
        return np.exp(-2j * np.pi * self.F_Hz * self.delay_s)

    @depB_property
    def h_a(self):
        self.dependencies('poles', 'a_vec')
        X = self.X_grid
        pval = self.get_default('poles', None)
        #pval = self.poles
        if pval is not None:
            val, lnG = pval.val_lnG(X)
            adj = self.phi_gain_adj_roots(pval)
            #print('Ap', np.max(abs(val)), lnG)
            return lnGxfer(adj * val * self._phasing_roots(pval), lnG)

        a_vec = self.get_default('a_vec', None)
        if a_vec is not None:
            X = self.X_grid
            val, lnG = self.poly.val_lnG(X, a_vec)
            adj = self.phi_gain_adj_pvec(a_vec)
            #print('Bp', np.max(abs(val)), lnG)
            return lnGxfer(adj * val * self._phasing_vec(a_vec), lnG)

        #otherwise compute it from the preferred poles representation
        pval = self.poles
        val, lnG = pval.val_lnG(X)
        adj = self.phi_gain_adj_roots(pval)
        #print('Cp', np.max(abs(val)), lnG)
        return lnGxfer(adj * val * self._phasing_roots(pval), lnG)

    @depB_property
    def h_b(self):
        self.dependencies('zeros', 'b_vec', 'lnG')
        X = self.X_grid
        zval = self.get_default('zeros', None)
        #zval = self.zeros
        if zval is not None:
            val, lnG = zval.val_lnG(X)
            adj = self.phi_gain_adj_roots(zval)
            #print('Az', np.max(abs(val)), lnG)
            return lnGxfer(adj * val * self._phasing_roots(zval), lnG)

        b_vec = self.get_default('b_vec', None)
        if b_vec is not None:
            val, lnG = self.poly.val_lnG(X, b_vec)
            adj = self.phi_gain_adj_pvec(b_vec)
            #print('Bz', np.max(abs(val)), lnG)
            return lnGxfer(adj * val * self._phasing_vec(b_vec), lnG)

        #otherwise compute it from the preferred poles representation
        zval = self.zeros
        val, lnG = zval.val_lnG(X)
        adj = self.phi_gain_adj_roots(zval)
        #print('Cz', np.max(abs(val)), lnG)
        return lnGxfer(adj * val * self._phasing_roots(zval), lnG)

    @depB_property
    def _extra_xfer_fit(self):
        if self.F_nyquist_Hz is None:
            xfer_p, lnG_p = self.poles_overlay.val_lnG(self.X_grid)
            xfer, lnG = self.zeros_overlay.val_lnG(self.X_grid, h = 1/xfer_p, lnG = -lnG_p)
        else:
            xfer_p, lnG_p = self.poles_overlay.val_lnG(self.Xn_grid)
            xfer, lnG = self.zeros_overlay.val_lnG(self.Xn_grid, h = 1/xfer_p, lnG = -lnG_p)
        return xfer * np.exp(lnG)

    @depB_property
    def data_effective(self):
        return self.data / (self._extra_xfer_fit * self.h_delay)

    @depB_property
    def xfer_fit(self):
        return self._extra_xfer_fit * (
            self.lnG.sign * np.exp(self.lnG.lnG - self.h_a.lnG + self.h_b.lnG)
            * (self.h_b.xfer * self.h_delay) / self.h_a.xfer
        )

    @depB_property
    def residuals(self):
        debias_reweight = 1/(.001 + self.W**2)
        retB = wavestate.bunch.Bunch()
        R = self.xfer_fit/self.data
        retB.resP = self.W * (R - 1)
        retB.resZ = self.W * (1/R - 1)
        retB.resD = self.W * (R - 1/R)
        retB.resD_average = np.sum(TFmath.abs_sq(retB.resD)) / (2 * len(self.data))
        retB.average = np.sum(
            (abs(retB.resP)**2 + abs(retB.resZ * debias_reweight)**2) / (1 + debias_reweight)
        ) / (2 * len(self.data))
        return retB

    @depB_property
    def residuals_average(self):
        return self.residuals.resD_average

    @depB_property
    def A_z(self):
        return self.V_b.vander * (
            (self._phasing_bfit() * self.lnG.sign) * np.exp(self.lnG.lnG - self.h_a.lnG + self.V_b.lnG)
            * self.W / (self.data_effective) / self.h_a.xfer
        ).reshape(-1, 1)

    @depB_property
    def A_zp(self):
        return self.V_a.vander * (
            (self._phasing_bfit() * self.lnG.sign) * np.exp(self.lnG.lnG - 2*self.h_a.lnG + self.V_a.lnG + self.h_b.lnG)
            * self.W / (self.data_effective) * self.h_b.xfer / self.h_a.xfer**2
        ).reshape(-1, 1)

    @depB_property
    def A_p(self):
        return self.V_a.vander * (
            (self._phasing_afit() * self.lnG.sign) * np.exp(-self.lnG.lnG - self.h_b.lnG + self.V_a.lnG)
            * self.W * (self.data_effective) / self.h_b.xfer
        ).reshape(-1, 1)

    @depB_property
    def A_pz(self):
        return self.V_b.vander * (
            (self._phasing_afit() * self.lnG.sign) * np.exp(-self.lnG.lnG - 2*self.h_b.lnG + self.h_a.lnG + self.V_b.lnG)
            * self.W * (self.data_effective) * self.h_a.xfer / self.h_b.xfer**2
        ).reshape(-1, 1)
