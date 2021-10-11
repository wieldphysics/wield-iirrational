# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
import scipy
#import declarative
import scipy.linalg
import scipy.special

from ..svd import SVD_SV
from .. import TFmath
from ..representations.polynomials import poly_constraints
from .rational_bases import DataFilterBase


def lstsq_safe(a, b, *args, **kwargs):
    """
    On some machines, LAPACK can fail on the gelsd driver call on some fits
    The other drivers are more resilient
    """
    driver = kwargs.get('lapack_driver', None)
    if driver is None:
        kwargs['lapack_driver'] = 'gelsd'
    try:
        return scipy.linalg.lstsq(a, b, *args, **kwargs)
    except Exception as E:
        import warnings
        warnings.warn("LAPACK lstsq error: {}\ntrying driver gelss".format(E))

    try:
        kwargs['lapack_driver'] = 'gelss'
        return scipy.linalg.lstsq(a, b, *args, **kwargs)
    except Exception as E:
        import warnings
        warnings.warn("LAPACK lstsq error: {}\ntrying driver gelsy".format(E))

    kwargs['lapack_driver'] = 'gelsy'
    return scipy.linalg.lstsq(a, b, *args, **kwargs)


class PolyFilterAlgorithmsBase(DataFilterBase):
    def match_pair_iter(
        self,
        Q_rank_cutoff,
        num_sequence = 2,
        zeros_first = True,
    ):
        npoles = self.npoles
        nzeros = self.nzeros

        self.matched_pairs_clear(Q_rank_cutoff = Q_rank_cutoff)

        self.npoles = len(self.poles)
        self.nzeros = len(self.zeros)

        if not zeros_first:
            self.fit_poles()
            num_sequence -= 1
        for idx in range(num_sequence):
            self.fit_zeros()
            self.fit_poles()

        if not zeros_first:
            self.fit_zeros()

        self.npoles = npoles
        self.nzeros = nzeros

    def matched_pairs_clear(
        self,
        Q_rank_cutoff = .5,
    ):
        """
        Match unique closest pairs, if they are within a bandwidth of 0Hz, then they are ignored
        """
        if Q_rank_cutoff == 0:
            return

        #cast into native representation
        poles = self.phi_rep_Snative(self.poles)
        zeros = self.phi_rep_Snative(self.zeros)

        #TODO, not working for rational_disc since it add branch-cut info
        #poles = self.RBalgo.expect(poles, self.RBalgo.root_constraints.mirror_real)
        #zeros = self.RBalgo.expect(zeros, self.RBalgo.root_constraints.mirror_real)

        rpB = TFmath.nearest_unique_pairs(zeros.c, poles.c)

        rpB.r12_list
        rpB.l1_remain
        rpB.l2_remain

        poles_new = []
        zeros_new = []
        for z, p in rpB.r12_list:
            Q_rank = self.Q_rank_Snative(z, p)
            if Q_rank > Q_rank_cutoff:
                poles_new.append(p)
                zeros_new.append(z)

        poles.c = np.concatenate([poles_new, rpB.l2_remain])
        zeros.c = np.concatenate([zeros_new, rpB.l1_remain])

        self.zeros = self.phi_Snative_rep(zeros)
        self.poles = self.phi_Snative_rep(poles)
        return

    def root_stabilize_Snative_MIP(
        self, rB,
        method = None,
    ):
        """
        MIP stands for modify-in-place, which must occur for the Snative
        representation
        """
        if method is None:
            return rB

        select_c = rB.c.real > 0
        select_r = rB.r.real > 0
        if method == 'flip':
            rB.r[select_r] = -rB.r[select_r]
            rB.c[select_c] = -rB.c[select_c].conjugate()
        elif method == 'remove':
            rB.r = rB.r[~select_r]
            rB.c = rB.c[~select_c]
        else:
            raise RuntimeError("Bad Argument")

        rB.clear()
        return

    def root_stabilize(self, rB, method = None):
        rB_native = self.phi_rep_Snative(rB)
        self.root_stabilize_Snative_MIP(rB_native, method = method)
        return self.phi_Snative_rep(rB_native)

    def Q_rank_Snative(self, p, z):
        z_soft = self.root_BW_soften_Snative(z)
        p_soft = self.root_BW_soften_Snative(p)
        return abs(p_soft - z_soft) * (1/p_soft.real**2 + 1/z_soft.real**2)**.5

    def root_BW_soften_Snative(self, root):
        #TODO, convert this into a distance_limit type function, to use the same
        #as ZPK_fitters
        F = root.imag
        F_abs = abs(F)
        BW = root.real
        BW_abs = abs(BW)

        if F_abs > self.F_max_Hz:
            dist = F_abs - self.F_max_Hz
            if BW_abs < dist/2:
                BW = np.sign(BW) * dist/2
        else:
            idx_nearest = np.searchsorted(self.F_Hz, F_abs)
            idx_min = max(idx_nearest - 1, 0)
            idx_max = min(idx_nearest + 2, len(self.F_Hz) - 1)
            BW_near = self.F_Hz[idx_max] - self.F_Hz[idx_min]
            if BW_abs < BW_near:
                BW = np.sign(BW) * BW_near
        return BW + 1j*F


class PolyFilterAlgorithms(PolyFilterAlgorithmsBase):
    """
    Mixin Class
    """
    poly_constraints = poly_constraints.eRoR

    def poly_fromroots(self, roots):
        """
        Cleans the coefficients as well
        """
        pvec, lnG = self.poly.fromroots_lnG(roots, X_scale = self.X_scale)
        return pvec.real

    def fit_poles_mod_zeros(
        self,
        **kwargs
    ):
        #remove the effect of linear (delay) phasing
        A_z = self.A_z
        #A_z = np.hstack([(1j*F_Hz).reshape(-1, 1), A_z])
        q, r = np.linalg.qr(A_z)

        A_p = self.A_p
        #the str is because py2 doesn't like some symbol in there
        A_p = A_p - np.einsum(str('ij,jk->ik'), q, np.einsum(str('ij,ik->jk'), q.conjugate(), A_p))

        S, V = SVD_SV(
            np.vstack([A_p.real, A_p.imag]),
            n_smallest  = 3,
            overwrite_a = True,
        )
        #print("POLES SVD: ", S[-4:])
        self.a_vec = V.T[:, -1]
        return

    def fit_SVD(
        self,
        modify = True,
        residuals = 0,
        **kwargs
    ):
        A_z = self.A_z
        A_p = self.A_p
        if residuals == 0:
            pass
        elif residuals > 0:
            A_z = self.V_b * (
                self.W  / (self.h_b / self.h_delay)
            ).reshape(-1, 1)
        elif residuals < 0:
            A_p = self.V_a * (
                self.W / (self.h_delay * self.h_a)
            ).reshape(-1, 1)

        A_svd = np.hstack([np.vstack([A_p.real, A_p.imag]), -np.vstack([A_z.real, A_z.imag])])

        S, V = SVD_SV(
            A_svd,
            n_smallest  = 3,
            overwrite_a = True,
        )

        #ignore division by zero warnings since NaN output is correct.
        #Argmin does the right thing with nan
        with np.errstate(divide='ignore', invalid='ignore'):
            S_rescaled = abs(S / V.T[0, :])

        idx_best = np.argmin(S_rescaled)
        if modify:
            self.a_vec = V.T[:A_p.shape[1], idx_best]
            self.b_vec = V.T[A_p.shape[1]:, idx_best]
        return S, V.T[0, :]

    def fit_zeros_mod_poles(
        self,
        **kwargs
    ):
        A_p = self.A_p
        #remove the effect of linear (delay) phasing
        #A_p = np.hstack([(1j*F_Hz).reshape(-1, 1), A_p])

        q, r = np.linalg.qr(A_p)

        A_z = self.A_z
        A_z = A_z - np.einsum(str('ij,jk->ik'), q, np.einsum(str('ij,ik->jk'), q.conjugate(), A_z))
        S, V = SVD_SV(
            np.vstack([A_z.real, A_z.imag]),
            n_smallest  = 3,
            overwrite_a = True,
        )
        #print("ZEROS SVD: ", S[-4:])
        self.b_vec = V.T[:, -1]
        return

    def fit_poles(self, **kwargs):
        #print(self.a_vec, self.b_vec)
        A_p = self.A_p
        #solve the problem with purely real taps
        a_fit, res, rank, s = lstsq_safe(
            np.vstack([A_p.real, A_p.imag]),
            np.hstack([self.W.real, self.W.imag]),
        )
        #print(a_fit)
        self.a_vec = a_fit
        return

    def fit_poles2(self, **kwargs):
        A_p = self.A_p
        #print(A_p.real.shape)
        a_fit, res, rank, s = lstsq_safe(
            np.block([[A_p.real, 0 * self.W.reshape(-1, 1)], [A_p.imag, (self.W * self.F_Hz).reshape(-1, 1)]]),
            np.hstack([self.W.real, self.W.imag]),
        )
        self.a_vec = a_fit[:-1]
        return

    def fit_zeros(self, **kwargs):
        A_z = self.A_z
        b_fit, res, rank, s = lstsq_safe(
            np.vstack([A_z.real, A_z.imag]),
            np.hstack([self.W.real, self.W.imag]),
        )
        self.b_vec = b_fit
        return

    def fit_zeros2(self, **kwargs):
        A_z = self.A_z
        b_fit, res, rank, s = lstsq_safe(
            np.block([[A_z.real, 0 * self.W.reshape(-1, 1)], [A_z.imag, (self.W * self.F_Hz).reshape(-1, 1)]]),
            np.hstack([self.W.real, self.W.imag]),
        )
        self.b_vec = b_fit[:-1]
        return


class PolyFilterAlgorithmsIm(PolyFilterAlgorithmsBase):
    """
    Mixin Class for polynomial algorithms which apply the constraint to mirror
    over the imaginary axis rather than the real.
    """
    poly_constraints = poly_constraints.eRoI

    def poly_fromroots(self, roots):
        """
        Cleans the coefficients as well
        """
        #print('ROOTS: ', roots)
        pvec, lnG = self.poly.fromroots_lnG(roots, X_scale = self.X_scale)
        if np.iscomplexobj(pvec):
            pvec[0::2].imag = 0
            pvec[1::2].real = 0
        return pvec

    #TODO
    def fit_SVD(
        self,
        modify = True,
        residuals = 0,
        **kwargs
    ):
        coeff_mult = np.empty(self.A_p.shape[1], dtype = complex)
        coeff_mult[0::2] = 1
        coeff_mult[1::2] = 1j
        A_p = coeff_mult.reshape(1, -1) * self.A_p

        coeff_mult = np.empty(self.A_z.shape[1], dtype = complex)
        coeff_mult[0::2] = 1
        coeff_mult[1::2] = 1j
        A_z = coeff_mult.reshape(1, -1) * self.A_z

        if residuals == 0:
            pass
        elif residuals > 0:
            A_z = self.V_b * (
                self.W  / (self.h_b / self.h_delay)
            ).reshape(-1, 1)
        elif residuals < 0:
            A_p = self.V_a * (
                self.W / (self.h_delay * self.h_a)
            ).reshape(-1, 1)

        A_svd = np.hstack([
            np.vstack([A_p.real, A_p.imag]),
            -np.vstack([A_z.real, A_z.imag]),
        ])

        S, V = SVD_SV(
            A_svd,
            n_smallest  = 3,
            overwrite_a = True,
        )
        with np.errstate(divide='ignore', invalid='ignore'):
            S_rescaled = abs(S / V.T[0, :])
        idx_best = np.argmin(S_rescaled)

        a_fit = V.T[:A_p.shape[1], idx_best]
        b_fit = V.T[A_p.shape[1]:, idx_best]

        a_fit = a_fit.astype(complex, copy = True)
        b_fit = b_fit.astype(complex, copy = True)

        a_fit[1::2] *= 1j
        b_fit[1::2] *= 1j

        if modify:
            self.a_vec = a_fit
            self.b_vec = b_fit
        return S, a_fit, b_fit

    #TODO
    def fit_poles(self, **kwargs):
        #print(self.a_vec, self.b_vec)
        coeff_mult = np.empty(self.A_p.shape[1], dtype = complex)
        coeff_mult[0::2] = 1
        coeff_mult[1::2] = 1j

        A_p = coeff_mult.reshape(1, -1) * self.A_p
        #solve the problem with purely real taps
        a_fit, res, rank, s = lstsq_safe(
            np.vstack([A_p.real, A_p.imag]),
            np.hstack([self.W.real, self.W.imag]),
        )
        #print("Sp: ", s / s[0])
        a_fit = a_fit.astype(complex, copy = True)
        a_fit[1::2] *= 1j
        self.a_vec = a_fit
        return

    #TODO
    def fit_zeros(self, **kwargs):
        coeff_mult = np.empty(self.A_z.shape[1], dtype = complex)
        coeff_mult[0::2] = 1
        coeff_mult[1::2] = 1j

        A_z = coeff_mult.reshape(1, -1) * self.A_z
        b_fit, res, rank, s = lstsq_safe(
            np.vstack([A_z.real, A_z.imag]),
            np.hstack([self.W.real, self.W.imag]),
        )
        #print("Sz: ", s / s[0])
        b_fit = b_fit.astype(complex, copy = True)
        b_fit[1::2] *= 1j
        self.b_vec = b_fit
        return

