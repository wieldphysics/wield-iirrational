# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
import scipy
import declarative
import scipy.linalg
import scipy.special

from .data_filtfit_split import DataFiltFitSplitBase
from ..svd import SVD_SV


def abs_sq(x):
    return x.real**2 + x.imag**2


class RationalDiscFilterRB(DataFiltFitSplitBase):
    phase_missing = False

    def __build__(
        self,
        _args = None,
        npoles      = None,
        nzeros      = None,
        ZPK         = ((), (), 1),
        mindelay    = False,
        stabilize   = False,
        parent      = None,
        poles_overlay = None,
        zeros_overlay = None,
        **kwargs
    ):
        if _args:
            raise RuntimeError("Only keyword arguments allowed")
        self.npoles    = npoles
        self.nzeros    = nzeros
        self.mindelay  = mindelay
        self.stabilize = stabilize

        if parent is not None:
            if poles_overlay is None:
                poles_overlay       = parent.poles_overlay
            if zeros_overlay is None:
                zeros_overlay       = parent.zeros_overlay
        else:
            if poles_overlay is None:
                poles_overlay = ()
            if zeros_overlay is None:
                zeros_overlay = ()

        self.poles_overlay = poles_overlay
        self.zeros_overlay = zeros_overlay

        super(RationalDiscFilterRB, self).__build__(
            parent = parent,
            **kwargs
        )

        #first to prevent lookups that shouldn't happen
        self.zeros     = ZPK[0]
        self.h_b_ref   = 1
        self.dependencies_for('zeros', 'b_vec')
        self.poles     = ZPK[1]
        self.h_a_ref   = 1
        self.dependencies_for('poles', 'a_vec')
        self.gain      = ZPK[2]

        @self.deco_generator(clear = False)
        def a_vec(self):
            """
            Setup so that the gain at the reference freq is ALWAYS 1
            """
            return self.poly.fromroots(self.poles).real

        @self.deco_setter(clear = False)
        def a_vec(self, val):
            gain = self.gain
            self.gain = gain / val[-1]
            self.h_a_ref = self.poly.val(self.Xp_ref_P, val) / val[-1]
            val = val / val[-1]
            if self.get_raw('stabilize'):
                poles = self.poly.roots(val)
                #use set raw so it won't update the gain redundantly/erroneously
                if np.any(abs(poles) > 1):
                    #only prevent using this vector if needed, stability is helped if we don't go convert between roots and vectors too much
                    #actually set it so that we get the stabilize action
                    self.poles = poles
                    return self.TAG_NO_SET
                else:
                    self.set_raw('poles', poles, setter = lambda s, x : x)

            self.dependencies('poles')
            return val

        @self.deco_generator(clear = False)
        def b_vec(self):
            """
            Setup so the gain is always specified at the ref
            """
            #TODO fix gain deps
            return self.get_raw('gain') * self.poly.fromroots(self.zeros).real

        @self.deco_setter(clear = False)
        def b_vec(self, val):
            val = np.asarray(val)
            gain = val[-1]
            self.gain = gain
            self.h_b_ref = self.poly.val(self.Xp_ref_Z, val) / val[-1]
            if self.get_raw('mindelay'):
                zeros = self.poly.roots(val)
                if np.any(abs(zeros) > 1):
                    #only prevent using this vector if needed, stability is helped if we don't go convert between roots and vectors too much
                    #actually set it so that we get the stabilize action
                    self.zeros = zeros
                    return self.TAG_NO_SET
                else:
                    #use set raw so it won't update the gain redundantly/erroneously
                    self.set_raw('zeros', zeros, setter = lambda s, x : x)
            #since they had to be computed anyway for the gain adjustment
            self.dependencies('zeros')
            return val

        @self.deco_generator(clear = False)
        def zeros(self):
            return self.poly.roots(self.b_vec)

        @self.deco_setter(clear = False)
        def zeros(self, val):
            #get previous
            val = np.asarray(val)
            if self.get_raw('mindelay'):
                val = self.root_stabilize(
                    val,
                    method = 'flip',
                    real_neg = 'ignore',
                    real_pos = 'ignore',
                )
            self.dependencies('b_vec')
            h_ref_now = self.poly.valfromroots(self.Xp_ref_Z, val)
            h_ref_prev = self.h_b_ref
            self.gain = self.get_raw('gain') * abs(h_ref_prev / h_ref_now)
            self.h_b_ref = h_ref_now
            return val

        @self.deco_generator(clear = False)
        def poles(self):
            return self.poly.roots(self.a_vec)

        @self.deco_setter(clear = False)
        def poles(self, val):
            val = np.asarray(val)
            if self.get_raw('stabilize'):
                val = self.root_stabilize(
                    val,
                    method = 'flip',
                    real_neg = 'ignore',
                    real_pos = 'ignore',
                )
            self.dependencies('a_vec')
            h_ref_now = self.poly.valfromroots(self.Xp_ref_P, val)
            h_ref_prev = self.h_a_ref
            self.gain = self.get_raw('gain') / abs(h_ref_prev / h_ref_now)
            self.h_a_ref = h_ref_now
            return val

        @self.deco_generator(clear = False)
        def zeros_full(self):
            return np.concatenate([self.zeros, self.zeros_overlay])

        @self.deco_generator(clear = False)
        def poles_full(self):
            return np.concatenate([self.poles, self.poles_overlay])

        @self.deco_setter(clear = False)
        def gain(self, val):
            #print("GAIN: ", val)
            return val

        @self.deco_generator(clear = False)
        def gain(self):
            raise RuntimeError("shouldn't compute")
            return self.b_vec[-1]/self.a_vec[-1]

        @self.deco_generator(clear = False)
        def V_b(self):
            return self.poly.vander(self.Xp_grid_Z, self.nzeros) * self.zeros_phasing.reshape(-1, 1)

        @self.deco_generator(clear = False)
        def V_a(self):
            return self.poly.vander(self.Xp_grid_P, self.npoles) * self.poles_phasing.reshape(-1, 1)

        @self.deco_generator(clear = False)
        def h_a(self):
            #TODO account for poles_phase_adjust
            self.dependencies('poles', 'a_vec')
            X = self.Xp_grid_P
            pval = self.get_default('poles', None)
            #pval = self.poles
            if pval is not None:
                val = self.poly.valfromroots(X, pval)
                return val * self.Xn_grid_P**len(pval)

            a_vec = self.get_default('a_vec', None)
            if a_vec is not None:
                X = self.Xp_grid_P
                val = self.poly.val(X, a_vec)
                return val * self.Xn_grid_P**(len(a_vec) - 1)

            #otherwise compute it from the preferred poles representation
            pval = self.poles
            val = self.poly.valfromroots(X, pval)
            return val * self.Xn_grid_P**len(pval)

        @self.deco_generator(clear = False)
        def h_b(self):
            #TODO account for zeros_phase_adjust
            self.dependencies('zeros', 'b_vec', 'gain')
            X = self.Xp_grid_Z
            zval = self.get_default('zeros', None)
            #zval = self.zeros
            if zval is not None:
                val = self.poly.valfromroots(X, zval)
                return self.gain * val * self.Xn_grid_Z**len(zval)

            b_vec = self.get_default('b_vec', None)
            if b_vec is not None:
                #don't compute and normalize by h_ref since the b-coeffs carry the gain
                val = self.poly.val(X, b_vec)
                return val * self.Xn_grid_Z**(len(b_vec) - 1)

            #otherwise compute it from the preferred poles representation
            zval = self.zeros
            val = self.poly.valfromroots(X, zval)
            return self.gain * val * self.Xn_grid_Z**len(zval)

        @self.deco_generator(clear = False)
        def _extra_xfer_fit(self):
            if not self.zeros_overlay and not self.poles_overlay:
                return 1
            return self.poly.valfromroots(self.Xp_grid_Z, self.zeros_overlay) / self.poly.valfromroots(self.Xp_grid_P, self.poles_overlay)

        @self.deco_generator(clear = False)
        def xfer_fit(self):
            return self._extra_xfer_fit * self.h_b / self.h_a

        @self.deco_generator(clear = False)
        def residuals(self):
            debias_reweight = 1/(.001 + self.W**2)
            retB = wavestate.bunch.Bunch()
            R = self.xfer_fit/self.data
            retB.resP = self.W * (R - 1)
            retB.resZ = self.W * (1/R - 1)
            retB.resD = self.W * (R - 1/R)
            retB.resD_average = np.sum(abs_sq(retB.resD)) / (2 * len(self.data))
            retB.average = np.sum(
                (abs(retB.resP)**2 + abs(retB.resZ * debias_reweight)**2) / (1 + debias_reweight)
            ) / (2 * len(self.data))
            return retB

        @self.deco_generator(clear = False)
        def residuals_average(self):
            return self.residuals.resD_average

        self.zeros_phase_adjust = None
        self.poles_phase_adjust = None

        @self.deco_generator(clear = False)
        def zeros_phasing(self):
            if self.zeros_phase_adjust is None:
                return self.Xn_grid_Z**self.nzeros
            elif self.zeros_phase_adjust < 0:
                return self.Xn_grid_Z**(self.nzeros + self.zeros_phase_adjust)
            else:
                return self.Xn_grid_Z**(self.zeros_phase_adjust)
            return

        @self.deco_generator(clear = False)
        def poles_phasing(self):
            if self.poles_phase_adjust is None:
                return self.Xn_grid_P**self.npoles
            elif self.poles_phase_adjust < 0:
                return self.Xn_grid_P**(self.npoles + self.poles_phase_adjust)
            else:
                return self.Xn_grid_P**(self.poles_phase_adjust)
            return

        @self.deco_generator(clear = False)
        def A_z(self):
            return self.V_b * (
                self.W / (self.data * self.h_a)
            ).reshape(-1, 1)

        @self.deco_generator(clear = False)
        def A_zp(self):
            return self.V_a * (
                self.W * self.h_b / (self.data * self.h_a**2)
            ).reshape(-1, 1)

        @self.deco_generator(clear = False)
        def A_p(self):
            return self.V_a * (
                self.W * (self.data / self.h_b)
            ).reshape(-1, 1)

        @self.deco_generator(clear = False)
        def A_pz(self):
            return self.V_b * (
                self.W * (self.data * self.h_a / self.h_b**2)
            ).reshape(-1, 1)

        return  # ~__init__

    @property
    def ZPK(self):
        #so cheap to make that we don't try to autodelete
        #TODO, fix gain term as it is currently based on reference frequency
        return (self.zeros, self.poles, self.gain)

    def root_stabilize(
        self,
        roots,
        method = None,
        real_neg = None,
        real_pos = None,
    ):
        roots = np.asarray(roots)
        select = (abs(roots) > 1)
        if real_neg is None:
            pass
        elif real_neg == 'ignore':
            select = select & ~((abs(roots.imag) < 1e-12) & (roots.real < 0))
        elif real_neg == 'ignore_float':
            select = select & ~((abs(roots.imag) < 1e-12) & (roots.real < 0))
            float_sel = ((abs(roots.imag) < 1e-12) & (roots.real < -1) & (roots.real > -5))
            roots[float_sel] = -5
        else:
            raise RuntimeError("Bad Argument")

        if real_pos is None:
            pass
        elif real_pos == 'ignore':
            select = select & ~((abs(roots.imag) < 1e-12) & (roots.real > 0))
        elif real_pos == 'ignore_float':
            select = select & ~((abs(roots.imag) < 1e-12) & (roots.real > 0))
            float_sel = ((abs(roots.imag) < 1e-12) & (roots.real > 1) & (roots.real < 5))
            roots[float_sel] = 5
        else:
            raise RuntimeError("Bad Argument")

        if method is None:
            return roots
        elif method == 'flip':
            roots = np.copy(roots)
            roots[select] = 1 / roots[select].conjugate()
            return roots
        elif method == 'remove':
            #just nullify them
            roots = roots[~select]
            return roots
        elif method == 'zero':
            roots[select] = 0
        else:
            raise RuntimeError("Bad Argument")
        return roots

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

    def matched_pairs_clear(
        self,
        Q_rank_cutoff = .5,
    ):
        """
        Match unique closest pairs, if they are within a bandwidth of 0Hz, then they are ignored
        """
        poles_r, poles_c  = self.poles_split
        zeros_r, zeros_c = self.zeros_split

        def nearest_idx(lst_1, lst_2):
            nearest_lst = []
            for r1 in lst_1:
                if r1 is None:
                    nearest_lst.append(None)
                    continue
                dist_nearest = float('inf')
                idx_nearest = None
                for idx_2, r2 in enumerate(lst_2):
                    if r2 is None:
                        continue
                    dist = abs(r1 - r2)
                    if dist < dist_nearest:
                        idx_nearest = idx_2
                        dist_nearest = dist
                nearest_lst.append(idx_nearest)
            return nearest_lst
        z_nearest = nearest_idx(zeros_c, poles_c)
        p_nearest = nearest_idx(poles_c, zeros_c)

        z_duals = []
        p_duals = []
        for idx_z, idx_p in enumerate(z_nearest):
            if idx_p is None:
                continue
            if idx_z != p_nearest[idx_p]:
                #not a unique pairing
                continue
            z = zeros_c[idx_z]
            p = poles_c[idx_p]
            Q_rank = abs(p-z) * (1/(1 - abs(p))**2 + 1/(1 - abs(z))**2)**.5
            if Q_rank < Q_rank_cutoff:
                z_duals.append(idx_z)
                p_duals.append(idx_p)
        p_duals = set(p_duals)
        z_duals = set(z_duals)

        poles_new = []
        for idx_p, pole in enumerate(poles_c):
            if idx_p in p_duals:
                continue
            poles_new.append(pole)
            poles_new.append(pole.conjugate())
        zeros_new = []
        for idx_z, zero in enumerate(zeros_c):
            if idx_z in z_duals:
                continue
            zeros_new.append(zero)
            zeros_new.append(zero.conjugate())
        poles_new.extend(poles_r)
        zeros_new.extend(zeros_r)
        self.poles = poles_new
        self.zeros = zeros_new
        return

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
        #remove the effect of linear (delay) phasing
        A_z = self.A_z
        A_p = self.A_p
        if residuals == 0:
            pass
        elif residuals > 0:
            A_z = self.V_b * (
                self.W * self.zeros_phasing / self.h_b
            ).reshape(-1, 1)
        elif residuals < 0:
            A_p = self.V_a * (
                self.W * self.poles_phasing / self.h_a
            ).reshape(-1, 1)

        A_svd = np.hstack([np.vstack([A_p.real, A_p.imag]), -np.vstack([A_z.real, A_z.imag])])

        S, V = SVD_SV(
            A_svd,
            n_smallest  = 3,
            overwrite_a = True,
        )
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
        a_fit, res, rank, s = scipy.linalg.lstsq(
            np.vstack([A_p.real, A_p.imag]),
            np.hstack([self.W.real, self.W.imag]),
        )
        #print(a_fit)
        print(a_fit.shape)
        self.a_vec = a_fit
        return

    def fit_zeros(self, **kwargs):
        A_z = self.A_z
        b_fit, res, rank, s = scipy.linalg.lstsq(
            np.vstack([A_z.real, A_z.imag]),
            np.hstack([self.W.real, self.W.imag]),
        )
        self.b_vec = b_fit
        return

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

    def fit_poles_s(self, **kwargs):
        A_p = self.A_p
        pt = np.exp(-1j * np.pi * np.linspace(self.F_Hz[-1], self.Fp_nyquist_Hz, 500))
        #pt = -1
        V_p = self.poly.vander(pt, self.npoles)
        prev = np.dot(V_p, list(self.a_vec) + [0] * (1 + self.nzeros - len(self.a_vec)))
        V_p = V_p / prev.reshape(-1, 1)
        W = 10
        a_fit, res, rank, s = scipy.linalg.lstsq(
            np.block([[A_p.real], [A_p.imag], [W * V_p.real], [W * V_p.imag]]),
            np.hstack([self.W.real, self.W.imag] + [W * np.ones(500)] + [np.zeros(500)]),
        )
        self.a_vec = a_fit[:]
        return


    def fit_poles2(self, **kwargs):
        #print(self.a_vec, self.b_vec)
        A_p = self.A_p
        #solve the problem with purely real taps
        print(A_p.real.shape)
        a_fit, res, rank, s = scipy.linalg.lstsq(
            #np.block([[A_p.real, 0 * self.W.reshape(-1, 1)], [A_p.imag, (self.W * 1j * self.F_Hz).reshape(-1, 1)]]),
            np.block([[A_p.real, 0 * self.W.reshape(-1, 1)], [A_p.imag, (self.W * self.F_Hz).reshape(-1, 1)]]),
            np.hstack([self.W.real, self.W.imag]),
        )
        print(a_fit.shape)
        self.a_vec = a_fit[:-1]
        return

    def fit_zeros2(self, **kwargs):
        A_z = self.A_z
        b_fit, res, rank, s = scipy.linalg.lstsq(
            np.block([[A_z.real, 0 * self.W.reshape(-1, 1)], [A_z.imag, (self.W * self.F_Hz).reshape(-1, 1)]]),
            np.hstack([self.W.real, self.W.imag]),
        )
        self.b_vec = b_fit[:-1]
        return

    def fit_zeros_s(self, **kwargs):
        A_z = self.A_z
        pt = np.exp(-1j * np.pi * np.linspace(self.F_Hz[-1], self.Fz_nyquist_Hz, 500))
        #pt = -1
        V_z = self.poly.vander(pt, self.nzeros)
        prev = np.dot(V_z, list(self.b_vec) + [0] * (1 + self.nzeros - len(self.b_vec)))
        V_z = V_z / prev.reshape(-1, 1)
        W = 10
        print(V_z.shape)
        b_fit, res, rank, s = scipy.linalg.lstsq(
            np.block([[A_z.real], [A_z.imag], [W * V_z.real], [W * V_z.imag]]),
            np.hstack([self.W.real, self.W.imag] + [W * np.ones(500)] + [np.zeros(500)]),
        )
        self.b_vec = b_fit[:]
        return

    def fit_zeros_s2(self, **kwargs):
        A_z = self.A_z

        Xp_grid = np.exp(1j * np.pi * self.F_Hz / np.max(self.F_Hz))
        V_b2 = self.poly.vander(Xp_grid, self.nzeros)
        z_phasing = Xp_grid**-self.nzeros

        A_z2 = V_b2 * (
            self.W * z_phasing / (self.data * self.h_a)
        ).reshape(-1, 1)

        q, r = np.linalg.qr(A_z2)

        A_z_proj = np.einsum(str('ij,jk->ik'), q, np.einsum(str('ij,ik->jk'), q.conjugate(), A_z))
        W_proj = np.einsum(str('ij,j->i'), q, np.einsum(str('ij,i->j'), q.conjugate(), self.W))

        b_fit, res, rank, s = scipy.linalg.lstsq(
            np.block([[A_z_proj.real], [A_z_proj.imag]]),
            np.hstack([W_proj.real, W_proj.imag]),
        )
        self.b_vec = b_fit
        return

    def fit_poles_s2(self, **kwargs):
        A_p = self.A_p

        Xp_grid = np.exp(1j * np.pi * self.F_Hz / np.max(self.F_Hz))
        V_a2 = self.poly.vander(Xp_grid, self.npoles)
        p_phasing = Xp_grid**-self.npoles

        A_p2 = V_a2 * (
            self.W * p_phasing * (self.data / self.h_b)
        ).reshape(-1, 1)

        q, r = np.linalg.qr(A_p2)

        A_p_proj = np.einsum(str('ij,jk->ik'), q, np.einsum(str('ij,ik->jk'), q.conjugate(), A_p))
        W_proj = np.einsum(str('ij,j->i'), q, np.einsum(str('ij,i->j'), q.conjugate(), self.W))

        a_fit, res, rank, s = scipy.linalg.lstsq(
            np.block([[A_p_proj.real], [A_p_proj.imag]]),
            np.hstack([W_proj.real, W_proj.imag]),
        )
        self.a_vec = a_fit
        return
