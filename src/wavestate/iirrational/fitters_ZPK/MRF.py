# -*- coding: utf-8 -*-
"""
"""


import numpy as np
import itertools
import contextlib
import declarative
from declarative import (
    NOARG,
    depB_property,
)


from ..external import scipy_optimize
from .. import TFmath
from .. import representations
from ..utilities import ensure_aid

from .data_filtfit_base import (
    DataFiltFitZ,
    DataFiltFitBase,
    DataFiltFitSf
)


class MultiReprFilterBase(DataFiltFitBase):
    #TODO add a "use_phase" option (defaulting to True)
    phase_missing = False
    rcond_cov = 1e-5
    root_constraint = representations.root_constraints.mirror_real
    RBalgo = representations.RBAlgorithms(constraint_standard = root_constraint)

    def __build__(
        self,
        _args                  = None,
        parent                 = None,
        gain_coding            = None,
        delay_coding           = None,
        num_codings            = None,
        den_codings            = None,
        parameters             = None,
        poles_overlay          = None,
        zeros_overlay          = None,
        residuals_type         = None,
        codings_ignore         = None,
        residuals_log_im_scale = None,
        distance_limit_auto    = None,
        distance_limit_scale   = None,
        delay_s_min            = None,
        delay_s_max            = None,
        coding_map             = None,
        h_infinity             = None,
        h_infinity_deweight    = None,
        max_BW_Hz              = None,
        copy                   = None,
        **kwargs
    ):
        if _args is not None:
            raise RuntimeError("Only keyword arguments allowed")
        if copy is not None:
            parent = copy
        super(MultiReprFilterBase, self).__build__(
            parent = parent,
            copy   = copy,
            **kwargs
        )
        if parent is not None and isinstance(parent, MultiReprFilterBase):
            if parameters is None:
                parameters             = parent.parameters
            if gain_coding is None:
                gain_coding            = parent.gain_coding
            if delay_coding is None:
                delay_coding           = parent.delay_coding
            #TODO, use the copier
            if num_codings is None:
                num_codings            = parent.num_codings
            #TODO, use the copier
            if den_codings is None:
                den_codings            = parent.den_codings
            if poles_overlay is None:
                poles_overlay          = parent.poles_overlay
            if zeros_overlay is None:
                zeros_overlay          = parent.zeros_overlay
            if residuals_type is None:
                residuals_type         = parent.residuals_type
            if codings_ignore is None:
                arg_codings_ignore     = parent.codings_ignore
            if residuals_log_im_scale is None:
                residuals_log_im_scale = parent.residuals_log_im_scale
            if distance_limit_auto is None:
                distance_limit_auto    = parent.distance_limit_auto
            if distance_limit_scale is None:
                distance_limit_scale = parent.distance_limit_scale
            if delay_s_min is None:
                delay_s_min = parent.delay_s_min
            if delay_s_max is None:
                delay_s_max = parent.delay_s_max
            if coding_map is None:
                coding_map = parent.coding_map
            if h_infinity is None:
                h_infinity = parent.h_infinity
            if h_infinity_deweight is None:
                h_infinity_deweight = parent.h_infinity_deweight
            if max_BW_Hz is None:
                max_BW_Hz = parent.max_BW_Hz
        else:
            if gain_coding is None:
                raise RuntimeError("Must specify gain_coding")
            if num_codings is None:
                raise RuntimeError("Must specify num_codings")
            if den_codings is None:
                raise RuntimeError("Must specify den_codings")
            if poles_overlay is None:
                poles_overlay = ()
            if zeros_overlay is None:
                zeros_overlay = ()
            if residuals_type is None:
                residuals_type = 'log'
            if codings_ignore is None:
                arg_codings_ignore = set()
            if residuals_log_im_scale is None:
                residuals_log_im_scale = 1
            if distance_limit_auto is None:
                distance_limit_auto = 0
            if distance_limit_scale is None:
                distance_limit_scale = 1
        arg_gain_coding       = gain_coding
        arg_delay_coding      = delay_coding
        arg_den_codings       = den_codings
        arg_num_codings       = num_codings

        self.max_BW_Hz = max_BW_Hz
        #must adjust this whenever something happens to the codings
        self.codings_revision  = 0

        self.coding_map          = coding_map
        self.h_infinity          = h_infinity
        self.h_infinity_deweight = h_infinity_deweight

        #can be 0: never, 1: on call, or 2: on all parameter update calls (during optimization)
        self.distance_limit_auto = distance_limit_auto
        self.distance_limit_scale = distance_limit_scale

        self.gain_coding            = arg_gain_coding
        self.num_codings            = arg_num_codings
        self.den_codings            = arg_den_codings

        if arg_delay_coding is not None:
            self.delay_coding           = arg_delay_coding
        if delay_s_min is not None:
            self.delay_s_min = delay_s_min
        if delay_s_max is not None:
            self.delay_s_max = delay_s_max

        self.poles_overlay          = poles_overlay
        self.zeros_overlay          = zeros_overlay
        self.residuals_type         = residuals_type
        self.residuals_log_im_scale = residuals_log_im_scale
        #assign after the setter is built
        self.codings_ignore = arg_codings_ignore
        self.auto_jacobian  = False

        #the parameters are mutable objects that need copy rules
        @self.deco_copier
        def gain_coding(self, val):
            return val.clone(self)

        @self.deco_copier
        def num_codings(self, val):
            return [c.clone(self) for c in val]

        @self.deco_copier
        def den_codings(self, val):
            return [c.clone(self) for c in val]

    @depB_property
    def h_infinity(self, val = NOARG):
        if val is NOARG or val is None:
            val = 0
        if callable(val):
            val = val(len(self.F_Hz))
        val = np.asarray(val)
        return val

    @depB_property
    def h_infinity_deweight(self, val = NOARG):
        if val is NOARG or val is None:
            val = 0.1

        val = float(val)
        return val

    @depB_property(autodeps = False)
    def residuals_log_im_scale(self, val = 1):
        #need to set without putting a dep on this value (or it will reset)
        self.codings_revision  = 1 + self.get_raw("codings_revision")
        return np.asarray(val)

    @depB_property
    def delay_s_min(self, val = NOARG):
        if val is NOARG:
            #default min delay is 22.5deg loss by the last frequency point
            #val = -(np.pi/8) / self.F_max_Hz
            val = 0
        else:
            #need to set without putting a dep on this value (or it will reset)
            self.codings_revision  = 1 + self.get_raw("codings_revision")
        return val

    @depB_property
    def delay_s_max(self, val = NOARG):
        if val is NOARG:
            #default max delay is 10deg loss by the last frequency point
            val = (np.pi/16) / self.F_max_Hz
        else:
            #need to set without putting a dep on this value (or it will reset)
            self.codings_revision  = 1 + self.get_raw("codings_revision")
        return val

    @depB_property
    def max_BW_Hz(self, val = NOARG):
        if val is NOARG or val is None:
            val = 2 * self.F_max_Hz
        return val

    @depB_property(autodeps = False)
    def distance_limit_auto(self, val):
        assert(val in [0, 1, 2])
        #need to set without putting a dep on this value (or it will reset)
        self.codings_revision  = 1 + self.get_raw("codings_revision")
        return val

    @depB_property
    def zeros_overlay(self, val = NOARG):
        if val is NOARG:
            return self.RBalgo.expect([])
        else:
            rB = self.RBalgo.expect(val)
            return rB

    @depB_property
    def poles_overlay(self, val = NOARG):
        if val is NOARG:
            return self.RBalgo.expect([])
        else:
            rB = self.RBalgo.expect(val)
            return rB

    @depB_property
    def _extra_xfer_fit(self):
        #TODO must add pole-zero rephasing for Z filters
        xfer_p, lnG_p = self.poles_overlay.val_lnG(self.Xex_grid)
        xfer, lnG = self.zeros_overlay.val_lnG(self.Xex_grid, h = 1/xfer_p, lnG = -lnG_p)
        return xfer * np.exp(lnG)

    @depB_property
    def data_effective(self):
        return self.data

    @depB_property(autodeps = False)
    def gain_coding(self, val = NOARG):
        if val is NOARG:
            raise RuntimeError("Currently can't generate a gain coding")
        else:
            #need to set without putting a dep on this value (or it will reset)
            self.codings_revision  = 1 + self.get_raw("codings_revision")
            return val.clone(self)

    @depB_property(autodeps = False)
    def delay_coding(self, val = NOARG):
        if val is NOARG:
            return self.coding_map.delay(self)
        else:
            return val.clone(self)

    @depB_property(autodeps = False)
    def num_codings(self, val):
        self.codings_revision = self.codings_revision + 1
        return [c.clone(self) for c in val]

    @depB_property(autodeps = False)
    def den_codings(self, val):
        self.codings_revision = self.codings_revision + 1
        return [c.clone(self) for c in val]

    @depB_property
    def root_codings_set(self):
        return set(self.num_codings + self.den_codings)

    @depB_property
    def codings_ignore(self, val):
        return val

    @depB_property
    def codings_ignore_mapped(self):
        self.codings_revision = self.codings_revision + 1
        val = self.codings_ignore
        if val is None:
            val = set()
        else:
            coding_id_map = {
                self.gain_coding.coding_id : self.gain_coding
            }
            if self.delay_coding is not None:
                coding_id_map[self.delay_coding.coding_id] = self.delay_coding
            for coding in self.num_codings:
                coding_id_map[coding.coding_id] = coding
            for coding in self.den_codings:
                coding_id_map[coding.coding_id] = coding

            codings_ignore_mapped = set()
            for coding in val:
                mcoding = coding_id_map.get(coding.coding_id, None)
                if mcoding is not None:
                    codings_ignore_mapped.add(mcoding)
            val = codings_ignore_mapped
        return val

    @depB_property
    def coding_lists(self):
        num_active = []
        num_ignore = []
        den_active = []
        den_ignore = []
        if self.gain_coding not in self.codings_ignore_mapped:
            num_active.append(self.gain_coding)
        else:
            num_ignore.append(self.gain_coding)

        if (
            self.delay_coding is not None
            and self.delay_coding not in self.codings_ignore_mapped
        ):
            num_active.append(self.delay_coding)
        else:
            num_ignore.append(self.delay_coding)

        for idx, coding in enumerate(self.num_codings):
            if coding not in self.codings_ignore_mapped:
                num_active.append(coding)
            else:
                num_ignore.append(coding)
        for idx, coding in enumerate(self.den_codings):
            if coding not in self.codings_ignore_mapped:
                den_active.append(coding)
            else:
                den_ignore.append(coding)
        return wavestate.bunch.Bunch(
            num_active = num_active,
            num_ignore = num_ignore,
            den_active = den_active,
            den_ignore = den_ignore,
        )

    @depB_property
    def xfer_fit(self, val = NOARG):
        if val is NOARG:
            #access just for the dependency
            #print("NEW_XFER")
            self.codings_revision
            if self.auto_jacobian:
                xfer, jac = self.transfer_jacobian_generate()
                self.residuals_jacobian = jac
                #TODO hack
                self.codings_revision
                return xfer
            else:
                xfer = self.transfer_generate()
                #TODO hack
                self.codings_revision
                return xfer
        else:
            #access just for the dependency
            self.codings_revision
            return val

    @depB_property
    def xfer_inactive(self):
        #access just for the dependency
        self.codings_revision
        return self.transfer_generate(generate_inactive = True)

    @depB_property
    def residuals_jacobian(self, val = NOARG):
        if val is NOARG:
            #access just for the dependency
            #print("NEW_JAC")
            self.codings_revision
            xfer, jac = self.transfer_jacobian_generate()
            #TODO can do this assignment, but need to set the dependencies with it
            self.xfer_fit = xfer
            return np.asarray(jac)
        else:
            #access just for the dependency
            self.codings_revision
            return val

    @depB_property
    def parameters(self, params = NOARG):
        if params is NOARG:
            params = []
            N_params = 0
            for coding in self.coding_lists.num_active:
                N_params += coding.N_parameters
                c_params = np.asarray(coding.reduce())
                #TODO, move these checks into the parameters
                if any(c_params.imag != 0):
                    raise RuntimeError("Coding not using real parameters {0}".format(coding))
                params.extend(c_params.real)
            for coding in self.coding_lists.den_active:
                N_params += coding.N_parameters
                c_params = np.asarray(coding.reduce())
                #TODO, move these checks into the parameters
                if any(c_params.imag != 0):
                    raise RuntimeError("Coding not using real parameters {0}".format(coding))
                params.extend(c_params.real)
            #add the dependency
            self.dependencies('codings_revision')
            return np.asarray(params, dtype = float)
        else:
            prev_params = self.get_raw('parameters', None)
            if (prev_params is not None) and np.all(prev_params == params):
                #print("SAME")
                return params

            #access for the dependency, and update the value to clear other dependent values
            self.codings_revision = self.codings_revision + 1

            N_current = 0
            for coding in self.coding_lists.num_active:
                num_next = coding.N_parameters
                coding.update(*params[N_current:N_current + num_next])
                N_current += num_next
            for coding in self.coding_lists.den_active:
                num_next = coding.N_parameters
                coding.update(*params[N_current:N_current + num_next])
                N_current += num_next
            assert(N_current == len(params))
            return params

    @depB_property
    def coding_maps(self):
        #TODO, use coding_lists, but have to handle gain_coding
        mapping = []
        locations_map = {}
        parameters_map = {}
        N_current = 0
        if self.gain_coding not in self.codings_ignore_mapped:
            mapping.extend(self.gain_coding.N_parameters * [self.gain_coding])
            num_next = self.gain_coding.N_parameters
            parameters_map[self.gain_coding] = N_current
            N_current += num_next
            locations_map[self.gain_coding] = 'gain_coding'

        if self.delay_coding is not None and self.delay_coding not in self.codings_ignore_mapped:
            mapping.extend(self.delay_coding.N_parameters * [self.delay_coding])
            num_next = self.delay_coding.N_parameters
            parameters_map[self.delay_coding] = N_current
            N_current += num_next
            locations_map[self.delay_coding] = 'delay_coding'

        num_idx = N_current
        for idx, coding in enumerate(self.num_codings):
            if coding in self.codings_ignore_mapped:
                continue
            num_next = coding.N_parameters
            parameters_map[coding] = N_current
            N_current += num_next
            mapping.extend(coding.N_parameters * [coding])
            locations_map[coding] = ('num_codings', idx)
        den_idx = N_current
        for idx, coding in enumerate(self.den_codings):
            if coding in self.codings_ignore_mapped:
                continue
            num_next = coding.N_parameters
            parameters_map[coding] = N_current
            N_current += num_next
            mapping.extend(coding.N_parameters * [coding])
            locations_map[coding] = ('den_codings', idx)
        #access for the dependency
        self.codings_revision
        return wavestate.bunch.Bunch(
            p2c     = mapping,
            c2loc   = locations_map,
            c2p     = parameters_map,
            num_idx = num_idx,
            den_idx = den_idx,
        )

    @depB_property
    def zeros(self, val = NOARG):
        if val is not NOARG:
            raise RuntimeError("Zeros not Directly Settable")
        #access just for the dependency
        self.codings_revision
        roots_c = []
        for coding in self.num_codings:
            c = coding.roots_c()
            assert(np.all(np.imag(c) >= 0))
            roots_c.extend(c)
        return representations.RootBunch(
            c = roots_c,
            r = [r for coding in self.num_codings for r in coding.roots_r()],
            constraint = representations.root_constraints.mirror_real
        )

    @depB_property
    def poles(self, val = NOARG):
        if val is not NOARG:
            raise RuntimeError("Poles not Directly Settable")
        #access just for the dependency
        self.codings_revision
        roots_c = []
        for coding in self.den_codings:
            c = coding.roots_c()
            assert(np.all(np.imag(c) >= 0))
            roots_c.extend(c)
        return representations.RootBunch(
            c = roots_c,
            r = [r for coding in self.den_codings for r in coding.roots_r()],
            constraint = representations.root_constraints.mirror_real
        )

    @depB_property
    def _gain_effect(self):
        #TODO, add gain properties to not need the raw parameter
        self.codings_revision
        gain_effect = 1
        for coding in self.num_codings:
            gain_effect *= coding.gain_effect
        for coding in self.den_codings:
            gain_effect /= coding.gain_effect
        return gain_effect

    @depB_property
    def gain(self, gain = NOARG):
        if gain is NOARG:
            #TODO, add gain properties to not need the raw parameter
            self.codings_revision
            return self.gain_coding.gain * self._gain_effect
        else:
            #TODO, add gain properties to not need the raw parameter
            gain_effect = self._gain_effect
            self.gain_coding.gain = gain / gain_effect
            self.codings_revision = self.codings_revision + 1
            return gain

    @depB_property
    def delay_s(self, val = NOARG):
        if val is NOARG:
            #TODO, add gain properties to not need the raw parameter
            self.codings_revision
            if self.delay_coding is None:
                return 0
            else:
                return self.delay_coding.delay_s
        else:
            self.codings_revision = self.codings_revision + 1
            self.delay_coding.delay_s = val
        return val

    @depB_property
    def residuals(self):
        self.codings_revision
        retB = wavestate.bunch.Bunch()
        R = self.xfer_fit/self.data
        retB.resP = self.W * (R - 1)
        retB.resZ = self.W * (1/R - 1)
        retB.resD = self.W * (R - 1/R)

        R_abs_sq = TFmath.abs_sq(R)
        R_abs = R_abs_sq**.5
        log_re = self.W * np.log(R_abs)
        log_im = self.W * R.imag / R_abs
        retB.resL = log_re + 1j * log_im * self.residuals_log_im_scale
        return retB

    @depB_property
    def residuals_preferred(self):
        self.codings_revision
        return self.residuals_compute(self.data, self.xfer_fit, self.W)

    @depB_property
    def residuals_average(self):
        if np.any(self.residuals_log_im_scale >= 1e-12):
            return np.sum(TFmath.abs_sq(self.residuals_preferred)) / (2 * len(self.data) - len(self.parameters))
        else:
            return np.sum(TFmath.abs_sq(self.residuals_preferred)) / (len(self.data) - len(self.parameters))

    @depB_property
    def residuals_max(self):
        #includes adjustment for parameter number. Not a real thing statistically...
        return np.max(TFmath.abs_sq(self.residuals_preferred))

    @depB_property
    def residuals_med(self):
        #includes adjustment for parameter number. Not a real thing statistically...
        return np.median(TFmath.abs_sq(self.residuals_preferred))

    @depB_property
    def fisher_matrix(self):
        #TODO, compute this using QR decomposition first
        #followed perhaps by SVD
        pjac = self.residuals_jacobian
        phess = np.dot(pjac, pjac.conjugate().T).real
        return phess

    @depB_property
    def fisher_matrix_rel(self):
        #TODO, compute this using QR decomposition first
        #followed perhaps by SVD
        mat = np.copy(self.fisher_matrix)
        d = mat.diagonal()**.5

        #recondition to prevent NANs
        d[d < 1e-15] = 1

        mat = mat / d.reshape(-1, 1)
        mat = mat / d.reshape(1, -1)
        return mat

    @depB_property
    def cov_matrix(self):
        #cov = np.linalg.pinv(self.fisher_matrix, rcond=1e-10)
        #uses the reconditioned fisher_matrix. This allows it to keep even the very badly behaved
        #values
        d = self.fisher_matrix.diagonal()**.5
        #recondition to prevent NANs
        d[d < 1e-15] = 1
        #TODO, find optimal rcond!
        p = np.linalg.pinv(
            self.fisher_matrix_rel,
            rcond = self.rcond_cov,
        )
        cov = p / d.reshape(1, -1) / d.reshape(-1, 1)
        #uses only the length of resP, so resZ would be the same
        return cov * self.residuals_average / (len(self.residuals_preferred) * 2 - cov.shape[0])

    @depB_property
    def cor_matrix(self):
        cov = np.copy(self.cov_matrix)
        d = cov.diagonal()**.5
        cov = cov / d.reshape(-1, 1)
        cov = cov / d.reshape(1, -1)
        return cov

    @depB_property
    def xfer_fit_error(self):
        cov = self.cov_matrix
        xfer, dfjac = self.transfer_jacobian_generate(xfer_fit_jacobian = True)
        dfjac = np.asarray(dfjac)
        perr_r = np.einsum(str("ij,ii,ij->j"), dfjac.real, cov, dfjac.real)
        perr_i = np.einsum(str("ij,ii,ij->j"), dfjac.imag, cov, dfjac.imag)

        dfjacAmp = np.einsum(str("ij,j->ij"), dfjac, abs(xfer) / xfer)
        perr_a = np.einsum(str("ij,ii,ij->j"), dfjacAmp.real, cov, dfjacAmp.real)
        dfjacRel = np.einsum(str("ij,j->ij"), dfjac, 1 / xfer)
        perr_ar = np.einsum(str("ij,ii,ij->j"), dfjacRel.real, cov, dfjacRel.real)
        perr_p = np.einsum(str("ij,ii,ij->j"), dfjacRel.imag, cov, dfjacRel.imag)
        cplx = perr_r**.5 + 1j * perr_i**.5
        return wavestate.bunch.Bunch(
            cplx  = cplx,
            mag   = perr_a**.5,
            mag_rel  = perr_ar**.5,
            mag_rel2 = perr_a**.5 / abs(xfer),
            phase = perr_p**.5,
        )

    @depB_property
    def residuals_error(self):
        cov = self.cov_matrix
        res = self.residuals_preferred
        dfjac = self.residuals_jacobian
        dfjac = np.asarray(dfjac)
        perr_r = np.einsum(str("ij,ii,ij->j"), dfjac.real, cov, dfjac.real)
        perr_i = np.einsum(str("ij,ii,ij->j"), dfjac.imag, cov, dfjac.imag)

        dfjacAmp = np.einsum(str("ij,j->ij"), dfjac, abs(res) / res)
        perr_a = np.einsum(str("ij,ii,ij->j"), dfjacAmp.real, cov, dfjacAmp.real)
        dfjacRel = np.einsum(str("ij,j->ij"), dfjac, 1 / res)
        perr_ar = np.einsum(str("ij,ii,ij->j"), dfjacRel.real, cov, dfjacRel.real)
        perr_p = np.einsum(str("ij,ii,ij->j"), dfjacRel.imag, cov, dfjacRel.imag)
        cplx = perr_r**.5 + 1j * perr_i**.5
        return wavestate.bunch.Bunch(
            cplx  = cplx,
            mag   = perr_a**.5,
            mag_rel   = perr_ar**.5,
            phase = perr_p**.5,
        )

    def transfer_generate(self, generate_inactive = False):
        if generate_inactive:
            #I believe these need a phasing adjustment if poles_overlay and zeros_overlay have different sizes
            #TODO also, if there are many poles_overlay and zeros_overlay, this is not numerically stable (should interlace)
            #TODO, should be using Xn_grid in the Z-domain
            xfer = self._extra_xfer_fit
            gain_effect = 1
            for coding_n, coding_d in itertools.zip_longest(
                    self.coding_lists.num_ignore,
                    self.coding_lists.den_ignore,
                    fillvalue = None,
            ):
                if coding_n is not None:
                    xfer = xfer * coding_n.transfer()
                    gain_effect *= coding_n.gain_effect
                if coding_d is not None:
                    xfer = xfer / coding_d.transfer()
                    gain_effect /= coding_d.gain_effect
            return xfer * gain_effect
        else:
            xfer = self.xfer_inactive
            gain_effect = 1
            for coding_n, coding_d in itertools.zip_longest(
                    self.coding_lists.num_active,
                    self.coding_lists.den_active,
                    fillvalue = None,
            ):
                if coding_n is not None:
                    xfer = xfer * coding_n.transfer()
                    gain_effect *= coding_n.gain_effect
                if coding_d is not None:
                    xfer = xfer / coding_d.transfer()
                    gain_effect /= coding_d.gain_effect
            return xfer * gain_effect

    def copy(self, **kwargs):
        return self.__class__(
            parent = self,
            **kwargs
        )

    def residuals_compute(
            self,
            data  = None,
            xfer  = None,
            W     = None,
            rtype = None
    ):
        if data is None:
            data = self.data
        if xfer is None:
            xfer = self.xfer_fit
        if W is None:
            W = self.W
        R = xfer/data
        return self.residuals_NLmap(R, W = W, rtype = rtype)

    def residuals_NLmap(self, R, W, rtype = None):
        if rtype is None:
            rtype = self.residuals_type
        if rtype == 'zeros':
            return W * (R - 1)
        elif rtype == 'poles':
            return W * (1/R - 1)
        elif rtype == 'dualA':
            return W * (R/2 + 0.5/R - 1)
        elif rtype == 'dualB':
            return W * (R - 1/R)/2
        elif rtype == 'log':
            R_abs_sq = TFmath.abs_sq(R)
            R_abs = R_abs_sq**.5
            select = (R_abs <= 0) | ~np.isfinite(R_abs)

            if any(select):
                #TODO, make this report better, requires logging or some kind
                #of hook
                #import sys
                #print(
                #    "BADNESS AT: ",
                #    self.F_Hz[select],
                #    R_abs[select],
                #    np.sum(select) / len(select),
                #    file = sys.stderr
                #)
                #print(
                #    self.gain_coding.gain,
                #)
                pass

            with np.errstate(divide='ignore', invalid='ignore'):
                log_re = W * np.log(R_abs)
                #TODO, test if there is a way to scale the residuals_log using the real part
                #log_resq = np.sum(log_re**2)
                log_im = W * R.imag / R_abs
            return log_re + 1j * log_im * self.residuals_log_im_scale
        return

    def transfer_jacobian_generate(
        self,
        xfer_fit_jacobian = False,
    ):
        #in principle it is faster to generate the xfer with the jacobian,
        #I am separating them as a wash
        # NEEDS SAME ORDERING AS PARAMETERS LIST

        #use separate bins for num and den so that the interlacing doesn't screw up the order
        op = wavestate.bunch.Bunch(
            jac_list_num = [],
            jac_list_den = [],
            xfer_running = 1
        )
        op.xfer_running = self.xfer_inactive
        rtype = self.residuals_type

        if xfer_fit_jacobian or rtype in ['poles', 'zeros', 'dualB', 'dualA']:
            def jac_modify(pjac, xfer, num):
                def generate(XY):
                    if num:
                        return XY * pjac / xfer
                    else:
                        return -XY * pjac / xfer
                return generate

        elif rtype == 'log':
            def jac_modify(pjac, xfer, num):
                jac_abs_sq = TFmath.abs_sq(xfer)
                jac_re = (pjac.real * xfer.real + pjac.imag * xfer.imag) / jac_abs_sq

                def generate(XY):
                    Y = XY / xfer
                    XY_abs = abs(XY)
                    jac_im = self.residuals_log_im_scale * ((pjac.imag * Y.real + pjac.real * Y.imag) - XY.imag * jac_re) / XY_abs
                    if num:
                        return jac_re + 1j * jac_im
                    else:
                        return -(jac_re + 1j * jac_im)
                return generate

        def coding_jac_append(coding, num):
            #don't include in jacobian if the codings_ignore contains the coding
            xfer, jac = coding.derivative_wtrans()

            if num:
                op.xfer_running = op.xfer_running * xfer * coding.gain_effect
                for pjac in jac:
                    op.jac_list_num.append(jac_modify(pjac, xfer, num = True))
            else:
                op.xfer_running = op.xfer_running / xfer / coding.gain_effect
                for pjac in jac:
                    op.jac_list_den.append(jac_modify(pjac, xfer, num = False))

        #interlace to prevent the value becoming crazy small and numerically unstable
        for coding_n, coding_d in itertools.zip_longest(
                self.coding_lists.num_active,
                self.coding_lists.den_active,
                fillvalue = None,
        ):
            if coding_n is not None:
                coding_jac_append(coding_n, num = True)
            if coding_d is not None:
                coding_jac_append(coding_d, num = False)

        xfer = op.xfer_running
        #have to do this after, otherwise they will be out of order from the interlacing
        jac_list = op.jac_list_num + op.jac_list_den
        #now need to actually inject the jacobian side of the chain rule, with the proper signs for the num, denom
        #first set are numerators of the xfer, including gain_coding
        if not xfer_fit_jacobian:
            if rtype == 'poles':
                XY = (xfer/self.data_effective)
            elif rtype == 'zeros':
                #minus sign to account for xfer on the denominator
                XY = -(self.data_effective/xfer)
            elif rtype == 'dualA':
                #minus sign goes to plus on the second term
                XY = (xfer/self.data_effective - self.data_effective/xfer)/2
            elif rtype == 'dualB':
                #minus sign goes to plus on the second term
                XY = (xfer/self.data_effective + self.data_effective/xfer)
            elif rtype == 'log':
                #minus sign goes to plus on the second term
                XY = xfer/self.data_effective
        else:
            XY = xfer

        for idx in range(len(jac_list)):
            jac_list[idx] = self.W * jac_list[idx](XY)

        return xfer, jac_list

    def optimize(
        self,
        **kwargs
    ):
        #start with Nelder-Mead as it is fast, not needing the jacobian, then use Levenberg-Marquardt to seal the deal
        #self.optimize_NM(**kwargs)
        return self.optimize_LS(**kwargs)

    opt_xtol_default    = 3e-16
    opt_ftol_default    = 3e-5
    opt_gtol_default    = 3e-5
    opt_method_default  = None
    opt_rescale_default = True
    opt_nfev_default    = 1000
    opt_log_default     = False

    def optimize_LS(
        self,
        xtol           = None,
        ftol           = None,
        gtol           = None,
        x_scale        = 'jac',
        max_nfev       = None,
        residuals_type = None,
        rescale        = True,
        aid            = None,
        log            = None,
        **kwargs
    ):
        aid = ensure_aid(aid)

        args = [
            'optimize_LS_{arg}.MRF',
            'optimize_{arg}.MRF',
            'optimize_LS_{arg}',
            'optimize_{arg}',
        ]
        xtol     = aid.hint_arg(xtol,     args, default = self.opt_xtol_default, arg = 'xtol')
        ftol     = aid.hint_arg(ftol,     args, default = self.opt_ftol_default, arg = 'ftol')
        gtol     = aid.hint_arg(gtol,     args, default = self.opt_gtol_default, arg = 'gtol')
        max_nfev = aid.hint_arg(max_nfev, args, default = self.opt_nfev_default, arg = 'max_nfev')

        rescale = rescale if rescale is not None else self.opt_rescale_default

        res_type_prev = self.residuals_type
        if residuals_type is not None:
            self.residuals_type = residuals_type
        if rescale:
            p_rescale = np.abs(np.asarray(self.parameters)) + 1e-5
        else:
            p_rescale = np.asarray(1)
        #TODO consider using auto_jacobian, gotta profile to see if it helps
        if np.any(abs(self.residuals_log_im_scale) > 1e-12):
            def opt(p):
                self.parameters = p * p_rescale
                return np.concatenate([
                    self.residuals_preferred.real,
                    self.residuals_preferred.imag,
                ])

            def opt_j(p, f = None):
                self.parameters = p * p_rescale
                jac = self.residuals_jacobian * p_rescale.reshape(-1, 1)
                return np.concatenate([
                    jac.real,
                    jac.imag
                ], axis = 1).T
        else:
            def opt(p):
                self.parameters = p * p_rescale
                return np.concatenate([
                    self.residuals_preferred.real,
                ])

            def opt_j(p, f = None):
                self.parameters = p * p_rescale
                jac = self.residuals_jacobian * p_rescale.reshape(-1, 1)
                return np.concatenate([
                    jac.real,
                ], axis = 1).T

        #if using h_infinity, then wrap the residuals again and apply the ranked mask
        if self.h_infinity.shape != () or self.h_infinity > 0:
            prev_opt = opt
            prev_opt_j = opt_j
            def opt(p):
                resid = prev_opt(p)
                resid_r = resid[:len(self.F_Hz)]
                resid_i = resid[len(self.F_Hz):]
                resid_sq = resid_r**2 + resid_i**2
                if self.h_infinity.shape == ():
                    if self.h_infinity_deweight == 0:
                        if self.h_infinity >= 1 - 1/len(self.F_Hz):
                            idx_max = np.argmax(resid_sq)
                            select = np.array([idx_max])
                        else:
                            argsort = np.argsort(-resid_sq)
                            select = argsort[:int((1 - self.h_infinity) * len(self.F_Hz) + 1)]
                        r = np.concatenate([resid_r[select], resid_i[select]])
                        return r
                    else:
                        argsort = np.argsort(-resid_sq)
                        select = argsort[:int((1 - self.h_infinity) * len(self.F_Hz) + 1)]
                        resid_r[select] *= self.h_infinity_deweight
                        resid_i[select] *= self.h_infinity_deweight
                        return np.concatenate([resid_r, resid_i])
                else:
                    argsort = np.argsort(-resid_sq)
                    resid_r[argsort] *= self.h_infinity
                    resid_i[argsort] *= self.h_infinity
                    return np.concatenate([resid_r, resid_i])

            def opt_j(p, f = None):
                resid    = prev_opt(p)
                resid_r  = resid[:len(self.F_Hz)]
                resid_i  = resid[len(self.F_Hz):]
                resid_sq = resid_r**2 + resid_i**2
                j = prev_opt_j(p, f)
                if self.h_infinity.shape == ():
                    if self.h_infinity_deweight == 0:
                        if self.h_infinity >= 1 - 1/len(self.F_Hz):
                            idx_max = np.argmax(resid_sq)
                            select = np.array([idx_max])
                        else:
                            argsort = np.argsort(-resid_sq)
                            select = argsort[:int((1 - self.h_infinity) * len(self.F_Hz) + 1)]
                        select = np.concatenate([select, len(self.F_Hz) + select])
                        return j[select, :]
                    else:
                        argsort = np.argsort(-resid_sq)
                        select = argsort[:int((1 - self.h_infinity) * len(self.F_Hz) + 1)]
                        j[select] *= self.h_infinity_deweight
                        j[len(self.F_Hz) + select] *= self.h_infinity_deweight
                        return j
                else:
                    argsort = np.argsort(-resid_sq)
                    j[argsort] *= self.h_infinity.reshape(-1, 1)
                    j[len(self.F_Hz) + argsort] *= self.h_infinity.reshape(-1, 1)
                    return j

        x0 = np.asarray(self.parameters) / p_rescale
        opt_r = scipy_optimize.trf(
            fun      = opt,
            jac      = opt_j,
            x0       = x0,
            J0       = opt_j(x0),
            f0       = opt(x0),
            xtol     = xtol,
            ftol     = ftol,
            gtol     = gtol,
            x_scale  = x_scale,
            max_nfev = max_nfev,
            verbose  = 0,
            **kwargs
        )

        #don't actually return an OptimizeResult because it kills matlab
        results = wavestate.bunch.Bunch()
        results.update(opt_r)
        self.parameters = results.x * p_rescale
        log = aid.hint_arg(log, args, default = self.opt_log_default, arg = 'log')
        if log:
            aid.log(results.message, results.nfev, results.njev)
        if residuals_type is not None:
            self.residuals_type = res_type_prev
        return results

    def optimize_NM(
        self,
        bootstrap      = True,
        tol            = 3e-16,
        xtol           = None,
        ftol           = None,
        method         = 'Nelder-Mead',
        max_nfev       = None,
        residuals_type = None,
        rescale        = True,
        aid            = None,
        **kwargs
    ):
        """
        Uses Nelder-Mead with a simplex that is bootstrapped from the
        svd of the jacobian. This gives particularly good
        initial directions to check
        """
        aid = ensure_aid(aid)
        if bootstrap:
            jac = np.concatenate([
                self.residuals_jacobian.real,
                self.residuals_jacobian.imag
            ], axis = 1)
            res = np.concatenate([
                self.residuals_preferred.real,
                self.residuals_preferred.imag
            ], axis = 0)
            d = np.sum(jac * jac, axis = 1)**.5
            #p_rescale the jacobian using the norms, makes the sum of the singular values equal to the number
            #of parameters (and greatly improves conditioning)
            d[abs(d) < 1e-10] = 1
            jac /= d.reshape(-1, 1)
            u, s, v = np.linalg.svd(
                jac,
                full_matrices=False
            )
            n_r = (res.dot(res))**.5
            v_r = v.dot(res) / n_r
            #deweight the smallest singular value as it will be replaced and it shouldn't screw up the fit if it is nasty
            #TODO, handle bad singular values
            u_proj = u * (v_r * n_r / (s + 1e-4 + np.min(s)) / d).reshape(-1, 1)
            u_best = np.sum(u_proj, axis = 1)
            #now replace the smallest singular vector with the optimal nearest quadratic solution
            u_proj[:, -1] = u_best
            initial_simplex = np.concatenate([np.zeros_like(d).reshape(-1, 1), u_proj], axis = 1).T + np.asarray(self.parameters)
        else:
            initial_simplex = None

        res_type_prev = self.residuals_type
        if residuals_type is not None:
            self.residuals_type = residuals_type
        try:

            args = [
                'optimize_NM_{arg}.MRF',
                'optimize_{arg}.MRF',
                'optimize_NM_{arg}',
                'optimize_{arg}',
            ]
            xtol     = aid.hint_arg(xtol,     args, default = self.opt_xtol_default, arg = 'xtol')
            ftol     = aid.hint_arg(ftol,     args, default = self.opt_ftol_default, arg = 'ftol')
            max_nfev = aid.hint_arg(max_nfev, args, default = self.opt_nfev_default, arg = 'max_nfev')
            rescale = rescale if rescale is not None else self.opt_rescale_default

            if rescale:
                p_rescale = (
                    np.abs(np.asarray(self.parameters))
                    + aid.hint('MRF_rescale_min', default = 1e-5)
                )
            else:
                p_rescale = np.asarray(1)

            if residuals_type is not None:
                self.residuals_type = residuals_type
            def opt(p):
                self.parameters = p * p_rescale
                return TFmath.norm_sq(self.residuals_preferred)

            if self.h_infinity.shape != () or self.h_infinity > 0:
                prev_opt = opt
                def opt(p):
                    resid = prev_opt(p)
                    resid_r = resid[:len(self.F_Hz)]
                    resid_i = resid[len(self.F_Hz):]
                    resid_sq = resid_r**2 + resid_i**2
                    if self.h_infinity.shape == ():
                        if self.h_infinity_deweight == 0:
                            if self.h_infinity >= 1 - 1/len(self.F_Hz):
                                idx_max = np.argmax(resid_sq)
                                select = np.array([idx_max])
                            else:
                                argsort = np.argsort(-resid_sq)
                                select = argsort[:int((1 - self.h_infinity) * len(self.F_Hz) + 1)]
                            r = np.concatenate([resid_r[select], resid_i[select]])
                            return r
                        else:
                            argsort = np.argsort(-resid_sq)
                            select = argsort[:int((1 - self.h_infinity) * len(self.F_Hz) + 1)]
                            resid_r[select] *= self.h_infinity_deweight
                            resid_i[select] *= self.h_infinity_deweight
                            return np.concatenate([resid_r, resid_i])
                    else:
                        argsort = np.argsort(-resid_sq)
                        resid_r[argsort] *= self.h_infinity
                        resid_i[argsort] *= self.h_infinity
                        return np.concatenate([resid_r, resid_i])

            kwargs['maxiter'] = max_nfev
            opt_r = scipy_optimize.neldermead(
                func            = opt,
                x0              = np.asarray(self.parameters) / p_rescale,
                xatol           = xtol,
                fatol           = ftol,
                maxiter         = max_nfev,
                maxfev          = max_nfev,
                initial_simplex = initial_simplex,
                callback        = None,
                disp            = False,
                return_all      = False,
            )
            results = wavestate.bunch.Bunch()
            results.update(opt_r)
            self.parameters = results.x * p_rescale

            if residuals_type is not None:
                self.residuals_type = res_type_prev

            return results
        finally:
            if residuals_type is not None:
                self.residuals_type = res_type_prev

    def optimize_min(
        self,
        tol            = 3e-16,
        xtol           = None,
        ftol           = None,
        method         = 'Nelder-Mead',
        max_nfev       = None,
        residuals_type = None,
        rescale        = True,
        aid            = None,
        **kwargs
    ):
        res_type_prev = self.residuals_type
        aid = ensure_aid(aid)

        args = [
            'optimize_min_{arg}.MRF',
            'optimize_{arg}.MRF',
            'optimize_min_{arg}',
            'optimize_{arg}',
        ]
        xtol     = aid.hint_arg(xtol,     args, default = self.opt_xtol_default, arg = 'xtol')
        ftol     = aid.hint_arg(ftol,     args, default = self.opt_ftol_default, arg = 'ftol')
        max_nfev = aid.hint_arg(max_nfev, args, default = self.opt_nfev_default, arg = 'max_nfev')

        method  = method if method is not None else self.opt_method_default
        rescale = rescale if rescale is not None else self.opt_rescale_default

        if rescale:
            p_rescale = (
                np.abs(np.asarray(self.parameters))
                + aid.hint('MRF_rescale_min', default = 1e-5)
            )
        else:
            p_rescale = np.asarray(1)

        if residuals_type is not None:
            self.residuals_type = residuals_type
        def opt(p):
            self.parameters = p * p_rescale
            return TFmath.norm_sq(self.residuals_preferred)

        def opt_j(p):
            self.parameters = p
            jac = self.residuals_jacobian * p_rescale.reshape(-1, 1)
            return 2 * (
                jac.real.dot(self.residuals_preferred.real)
                + jac.imag.dot(self.residuals_preferred.imag)
            )

        def opt_h(p):
            self.parameters = p
            jac = self.residuals_jacobian * p_rescale.reshape(-1, 1)
            hess = 2 * (
                jac.real.dot(jac.T.real)
                + jac.imag.dot(jac.T.imag)
            )
            return hess

        more_args = dict()
        if method in ['Nelder-Mead']:
            #doesn't need jac or hess
            pass
        else:
            more_args['jac'] = opt_j
            more_args['hess'] = opt_h

        kwargs['maxiter'] = max_nfev
        import scipy.optimize
        opt_r = scipy.optimize.minimize(
            opt,
            np.asarray(self.parameters) / p_rescale,
            tol     = tol,
            method  = method,
            options = kwargs,
            **more_args
        )
        results = wavestate.bunch.Bunch()
        results.update(opt_r)
        self.parameters = results.x * p_rescale

        if residuals_type is not None:
            self.residuals_type = res_type_prev

        return results

    def xfer_eval(self, F_Hz):
        return self.ZPKrep.xfer_eval(F_Hz = F_Hz,)

    def distance_limit(self, F_Hz, with_derivative = False):
        F_idx = np.searchsorted(self.F_Hz, F_Hz) - 1
        if F_idx < 1:
            val = self.F_Hz[1] - self.F_Hz[0]
            assert(val > 0)
            if not with_derivative:
                return val
            else:
                return val, 0
        if F_idx >= len(self.F_Hz) - 2:
            val = F_Hz - self.F_Hz[-4]
            #TODO
            #val = self.F_Hz[-1] - self.F_Hz[-4]
            assert(val > 0)
            if not with_derivative:
                return val
            else:
                return val, 1
        #get the local mean width of bins to restrict bandwidths
        F_mean_Hz_2 = (self.F_Hz[F_idx+2] - self.F_Hz[F_idx]) / 2
        F_mean_Hz_1 = (self.F_Hz[F_idx+1] - self.F_Hz[F_idx-1]) / 2
        x = (F_Hz - self.F_Hz[F_idx] ) / (self.F_Hz[F_idx + 1] - self.F_Hz[F_idx])
        if self.data is None:
            ratio = 1
        else:
            W = self.W
            if W is None:
                W = 1
            if np.asarray(W).shape != self.data.shape:
                W = 5
            else:
                W = W[F_idx-1:F_idx+3]
            dabs = abs(self.data[F_idx-1:F_idx+3])
            y_max = max(abs((1 - 1/W)) * dabs)
            y_min = min(abs((1 + 1/W)) * dabs)
            ratio = y_min/y_max
            if ratio > 1:
                ratio = 1
        #TODO
        #use the SNR to make a local density estimate!
        val = F_mean_Hz_1*(1 - x) + x * F_mean_Hz_2
        if (val < 0):
            print("VAL: ", F_mean_Hz_1, F_mean_Hz_2, x)
            assert(False)
        assert(self.distance_limit_scale > 0)
        assert(ratio > 0)

        if not with_derivative:
            return ratio * val * self.distance_limit_scale
        else:
            d = (F_mean_Hz_2 - F_mean_Hz_1) / (self.F_Hz[F_idx + 1] - self.F_Hz[F_idx])
            return ratio * val * self.distance_limit_scale, ratio * d * self.distance_limit_scale

    def optimize_delay(
        self,
        only_delay = True,
        **kwargs
    ):
        if self.delay_coding is None:
            return
        #store the previous so they may be restored
        prev_codings_ignore = self.codings_ignore

        if only_delay:
            #ignore all other codings
            self.codings_ignore = set(self.num_codings + self.den_codings)

        ret = self.optimize(**kwargs)

        #restore
        self.codings_ignore = prev_codings_ignore
        return ret

    def codings_in_region(
        self,
        F_min_Hz,
        F_max_Hz,
        num = True,
        den = True,
    ):
        raise NotImplementedError()
        #TODO, make compatible with S domain
        #TODO: add more tests and bandwidth limits
        gotten = []
        def test_get(coding):
            r_lst = coding.roots_c_Sf(self)
            if len(r_lst) == 1:
                r = r_lst[0]
                F = np.angle(r) / np.pi * self.F_nyquist_Hz
                if F_max_Hz is not None and F > F_max_Hz:
                    return
                if F_min_Hz is not None and F < F_min_Hz:
                    return
                gotten.append(coding)

        if num:
            for coding in self.num_codings:
                test_get(coding)
        if den:
            for coding in self.den_codings:
                test_get(coding)
        return gotten

    def regenerate(
        self,
        ZPKrep = None,
        zeros = None,
        poles = None,
        gain = None,
        p_c = None,
        z_c = None,
        p_r = None,
        z_r = None,
        coding_map = None,
        **kwargs
    ):
        from .ZPKrep2MRF import MRF2MRF
        return MRF2MRF(
            self,
            ZPKrep = ZPKrep,
            zeros  = zeros,
            poles  = poles,
            gain   = gain,
            p_c = p_c,
            z_c = z_c,
            p_r = p_r,
            z_r = z_r,
            coding_map = coding_map,
            **kwargs
        )

    @contextlib.contextmanager
    def with_codings_only(self, codings, prev_ignores = True):
        prev_ignores = self.codings_ignore
        if prev_ignores:
            self.codings_ignore = self.codings_ignore.union(self.root_codings_set - set(codings))
        else:
            self.codings_ignore = self.root_codings_set - set(codings)
        yield
        self.codings_ignore = prev_ignores
        return

    @contextlib.contextmanager
    def without_codings_only(self, codings, prev_ignores = True):
        prev_ignores = self.codings_ignore
        if prev_ignores:
            self.codings_ignore = self.codings_ignore.union(set(codings))
        else:
            self.codings_ignore = set(codings)
        yield
        self.codings_ignore = prev_ignores


def check_jac(fitterMRF, idx, scaling = .01):
    p_orig = np.asarray(fitterMRF.parameters)
    jac = fitterMRF.residuals_jacobian
    res1 = fitterMRF.residuals_preferred
    shift = [0] * len(p_orig)
    offset = scaling / (np.dot(jac[idx], jac[idx].conjugate()).real)**.5
    shift[idx] = offset
    fitterMRF.parameters = p_orig + np.asarray(shift)
    res2 = fitterMRF.residuals_preferred
    njac = (res2 - res1) / offset

    fitterMRF.parameters = p_orig
    with np.errstate(divide='ignore', invalid='ignore'):
        return njac / jac[idx], TFmath.abs_sq(jac[idx]) > 0


class MultiReprFilterZ(MultiReprFilterBase, DataFiltFitZ):
    pass


class MultiReprFilterS(MultiReprFilterBase, DataFiltFitSf):
    pass
