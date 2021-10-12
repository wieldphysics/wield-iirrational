# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
import scipy
import declarative
import scipy.linalg

from .data_filtfit_base import DataFiltFitBase
from .roots_bin import (roots_bin_palindromicz, roots_re_pair)
from .svd import SVD_SV


def abs_sq(x):
    return x.real**2 + x.imag**2


class RationalDiscFilterMag(DataFiltFitBase):
    phase_missing = False

    def __build__(
        self,
        _args  = None,
        npoles = None,
        nzeros = None,
        ZPK    = ((), (), 1),
        parent = None,
        **kwargs
    ):
        if _args:
            raise RuntimeError("Only keyword arguments allowed")
        self.npoles    = npoles
        self.nzeros    = nzeros

        super(RationalDiscFilterMag, self).__build__(
            parent = parent,
            **kwargs
        )

        #first to prevent lookups that shouldn't happen
        self.hsq_b_ref   = 1
        self.dependencies_for('zeros', 'bsq_vec')
        self.hsq_a_ref   = 1
        self.dependencies_for('poles', 'asq_vec')
        self.zeros     = ZPK[0]
        self.poles     = ZPK[1]

        self.gain_sq   = 1

        @self.deco_generator()
        def data_magsq(self):
            return abs_sq(self.data)

        @self.deco_generator(clear = False)
        def asq_vec(self):
            """
            Setup so that the gain at the reference freq is ALWAYS 1
            """
            vec = self.poly.fromroots(self.poles).real
            asq_vec_full = np.convolve(vec, vec[::-1] / vec[0], mode='full')
            asq_vec = asq_vec_full[len(vec)-1:]
            self.dependencies('poles')
            return asq_vec

        @self.deco_setter(clear = False)
        def asq_vec(self, val):
            gain_sq = self.gain_sq
            self.gain_sq = gain_sq / abs(val[-1])
            val = val / abs(val[-1])
            self.dependencies('poles')
            return val

        @self.deco_generator(clear = False)
        def bsq_vec(self):
            """
            Setup so the gain is always specified at the ref
            """
            #TODO fix gain deps
            vec = self.get_raw('gain') * self.poly.fromroots(self.zeros).real
            bsq_vec_full = np.convolve(vec, vec[::-1] / vec[0], mode='full')
            bsq_vec = bsq_vec_full[len(vec)-1:]
            self.dependencies('zeros')
            return bsq_vec

        @self.deco_setter(clear = False)
        def bsq_vec(self, val):
            val = np.asarray(val)
            self.gain_sq = abs(val[-1])
            val = val / abs(val[-1])
            self.dependencies('zeros')
            return val

        @self.deco_generator(clear = False)
        def zeros_split(self):
            vec = np.concatenate([self.bsq_vec[1::][::-1], self.bsq_vec])
            return roots_bin_palindromicz(self.poly.roots(vec), F_nyquist_Hz = self.F_nyquist_Hz)

        @self.deco_generator(clear = False)
        def zeros(self):
            self.dependencies('bsq_vec')
            return roots_re_pair(*self.zeros_split)

        @self.deco_setter(clear = False)
        def zeros(self, val):
            #get previous
            val = np.asarray(val)
            self.dependencies('bsq_vec')
            self.hsq_b_ref = self.poly.valfromroots(self.Xp_ref, val)
            return val

        @self.deco_generator
        def h_b_ref(self):
            return self.hsq_b_ref**.5

        @self.deco_generator
        def h_a_ref(self):
            return self.hsq_a_ref**.5

        @self.deco_generator(clear = False)
        def poles_split(self):
            vec = np.concatenate([self.asq_vec[1::][::-1], self.asq_vec])
            r_r, r_c = roots_bin_palindromicz(
                self.poly.roots(vec),
                F_nyquist_Hz = self.F_nyquist_Hz
            )
            return r_r, r_c

        @self.deco_generator(clear = False)
        def poles(self):
            self.dependencies('asq_vec')
            return roots_re_pair(*self.poles_split)

        @self.deco_setter(clear = False)
        def poles(self, val):
            val = np.asarray(val)
            self.dependencies('asq_vec')
            h_ref_now = abs_sq(self.poly.valfromroots(self.Xp_ref, val))
            h_ref_prev = self.hsq_a_ref
            #self.gain = self.get_raw('gain') / abs(h_ref_prev / h_ref_now)**.5
            self.hsq_a_ref = h_ref_now
            return val

        @self.deco_generator
        def poles_full(self):
            return np.asarray(self.poles)

        @self.deco_generator
        def zeros_full(self):
            return np.asarray(self.zeros)

        @self.deco_setter(clear = True)
        def gain(self, val):
            self.dependencies('bsq_vec')
            return val

        @self.deco_generator(clear = True)
        def gain(self):
            self.dependencies('bsq_vec')
            #TODO, could make more stable
            gain2 = (np.product(np.abs(self.poles)) / np.product(np.abs(self.zeros)))**.5
            return self.gain_sq**.5 * gain2
        self.gain      = ZPK[2]

        @self.deco_generator(clear = False)
        def V_b(self):
            #return self.poly.vander(self.Xc_grid, self.nzeros).real
            #TODO there is a more efficient way to construct this
            v = self.poly.vander(self.Xp_grid, self.nzeros).real
            v[:, 1:] *= 2
            return v

        @self.deco_generator(clear = False)
        def V_a(self):
            #return self.poly.vander(self.Xc_grid, self.npoles).real
            #TODO there is a more efficient way to construct this
            v = self.poly.vander(self.Xp_grid, self.npoles).real
            v[:, 1:] *= 2
            return v

        @self.deco_generator
        def V_ref_b(self):
            v = self.poly.vander(self.Xp_ref, self.nzeros).real
            v[:, 1:] *= 2
            return v

        @self.deco_generator
        def V_ref_a(self):
            v = self.poly.vander(self.Xp_ref, self.npoles).real
            v[:, 1:] *= 2
            return v

        @self.deco_generator(clear = False)
        def h_a(self):
            #TODO account for phase adjust
            self.dependencies('poles', 'asq_vec')
            val = self.poly.valfromroots(self.Xp_grid, self.poles)
            return val

        @self.deco_generator(clear = False)
        def h_b(self):
            #TODO account for phase adjust
            self.dependencies('zeros', 'bsq_vec', 'gain')
            val = self.poly.valfromroots(self.Xp_grid, self.zeros)
            return self.gain * val

        @self.deco_generator(clear = False)
        def hsq_a(self):
            self.dependencies('poles', 'asq_vec')
            return (np.dot(self.V_a[:, :len(self.asq_vec)], self.asq_vec))

        @self.deco_generator(clear = False)
        def hsq_b(self):
            #TODO account for zeros_phase_adjust
            self.dependencies('zeros', 'bsq_vec', 'gain')
            return (np.dot(self.V_b[:, :len(self.bsq_vec)], self.bsq_vec))

        @self.deco_generator(clear = False)
        def xfer_fit(self):
            return self.xfer_fit_magsq**.5
            return self.h_b / self.h_a

        @self.deco_generator(clear = False)
        def xfer_fit_magsq(self):
            return self.gain_sq * self.hsq_b / self.hsq_a

        @self.deco_generator(clear = False)
        def residuals(self):
            debias_reweight = 1/(.001 + self.W**2)
            retB = wavestate.bunch.Bunch()
            R = self.xfer_fit_magsq / self.data_magsq
            retB.resP = self.W * (R - 1) / 2
            retB.resZ = self.W * (1/R - 1) / 2
            retB.resD = self.W * (R - 1/R) / 2
            retB.resD_average = np.sum(abs(retB.resD)**2) / (4 * len(self.data))
            retB.average = np.sum(
                (abs(retB.resP)**2 + abs(retB.resZ * debias_reweight)**2) / (1 + debias_reweight)
            ) / (2 * len(self.data))
            return retB

        @self.deco_generator(clear = False)
        def residuals_average(self):
            return self.residuals.resD_average

        @self.deco_generator(clear = False)
        def A_z(self):
            return self.V_b * (
                self.W / (self.data_magsq  * self.hsq_a)
            ).reshape(-1, 1)

        @self.deco_generator(clear = False)
        def A_zp(self):
            return self.V_a * (
                self.W * self.hsq_b / (self.data_magsq * self.hsq_a**2)
            ).reshape(-1, 1)

        @self.deco_generator(clear = False)
        def A_p(self):
            return self.V_a * (
                self.W * (self.data_magsq / self.hsq_b)
            ).reshape(-1, 1)

        @self.deco_generator(clear = False)
        def A_pz(self):
            return self.V_b * (
                self.W * (self.data_magsq * self.hsq_a / self.hsq_b**2)
            ).reshape(-1, 1)
        return  # ~__init__

    def matched_pairs_clear(
        self,
        Q_rank_cutoff = .5,
    ):
        """
        Match unique closest pairs, if they are within a bandwidth of 0Hz, then they are ignored
        """
        poles_r, poles_c = self.poles_split
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
        max_size  = None,
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
            A_p.real,
            n_smallest  = 3,
            overwrite_a = True,
        )
        #print("POLES SVD: ", S[-4:])
        self.asq_vec = V.T[:, -1]
        return

    def fit_zeros_mod_poles(
        self,
        max_size  = None,
        **kwargs
    ):
        A_p = self.A_p
        #remove the effect of linear (delay) phasing
        #A_p = np.hstack([(1j*F_Hz).reshape(-1, 1), A_p])

        q, r = np.linalg.qr(A_p)

        A_z = self.A_z
        A_z = A_z - np.einsum(str('ij,jk->ik'), q, np.einsum(str('ij,ik->jk'), q.conjugate(), A_z))
        S, V = SVD_SV(
            A_z.real,
            n_smallest  = 3,
            overwrite_a = True,
        )
        #print("ZEROS SVD: ", S[-4:])
        self.bsq_vec = V.T[:, -1]
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

    def fit_poles(self, **kwargs):
        #print(self.asq_vec, self.bsq_vec)
        A_p = self.A_p
        #solve the problem with purely real taps
        a_fit, res, rank, s = scipy.linalg.lstsq(
            A_p.real,
            self.W.real,
        )
        self.asq_vec = a_fit
        return

    def fit_zeros(self, **kwargs):
        A_z = self.A_z
        b_fit, res, rank, s = scipy.linalg.lstsq(
            A_z.real,
            self.W.real,
        )
        self.bsq_vec = b_fit
        return

    def fit_polesX(self, **kwargs):
        #print(self.asq_vec, self.bsq_vec)
        A_p = self.A_p
        #solve the problem with purely real taps
        a_fit, res, rank, s = scipy.linalg.lstsq(
            self.data_magsq.reshape(-1, 1) * A_p.real,
            self.data_magsq * self.W.real,
        )
        self.asq_vec = a_fit
        return

    def fit_zerosX(self, **kwargs):
        A_z = self.A_z
        b_fit, res, rank, s = scipy.linalg.lstsq(
            1/self.data_magsq.reshape(-1, 1) * A_z.real,
            1/self.data_magsq * self.W.real,
        )
        self.bsq_vec = b_fit
        return

    def remove_doublets(self):
        p_r, p_c = self.poles_split
        z_r, z_c = self.zeros_split

        thresh_var = 1e-2
        thresh_covar = .1

        #mapping from zeros to poles, inside bunch with covar data
        used_zs = dict()

        for idx_p, p in enumerate(p_c):
            min_idx_z = None
            min_dist = 1
            for idx_z, z in enumerate(z_c):
                dist = abs_sq(p - z)
                if dist < min_dist:
                    min_idx_z = idx_z
                    min_dist = dist

            if min_idx_z is None:
                continue
            rp, rz, rcov_n = self.generate_ZPK_covar(
                cpoles_seq = [idx_p],
                czeros_seq = [min_idx_z]
            )
            if rp > thresh_var and rz > thresh_var and rcov_n > thresh_covar:
                if min_idx_z in used_zs:
                    if rcov_n < used_zs[min_idx_z].rcov_n:
                        #should maybe find second closest
                        continue

                #print(rp, rz, rcov_n, min_dist**.5)
                used_zs[min_idx_z] = wavestate.bunch.Bunch(
                    idx_p  = idx_p,
                    reff   = (rp**-2 + rz**-2)**(-.5),
                    rp     = rp,
                    rz     = rz,
                    rcov_n = rcov_n,
                )
        drop_ps = []
        drop_zs = []
        for idx_z, pB in used_zs.items():
            drop_zs.append(idx_z)
            drop_ps.append(pB.idx_p)

        drop_ps.sort()
        drop_zs.sort()

        def drop_join(roots, drop):
            c_r_join = []
            idx_r_prev = -1
            for idx_r in drop:
                c_r_join.append(roots[idx_r_prev+1:idx_r])
                idx_r_prev = idx_r
            c_r_join.append(roots[idx_r_prev+1:])
            from itertools import chain
            return tuple(chain(*c_r_join))

        p_c_new = drop_join(p_c, drop_ps)
        z_c_new = drop_join(z_c, drop_zs)
        self.poles = tuple(p_r) + p_c_new + tuple(r.conjugate() for r in p_c_new)
        self.zeros = tuple(z_r) + z_c_new + tuple(r.conjugate() for r in z_c_new)
        self.npoles = len(self.poles)
        self.nzeros = len(self.zeros)
        return

    def fit_pzpz(
        self,
        max_size          = None,
        collect_all       = False,
        zeros_first       = False,
        n_svd             = 1,
        n_iter            = 10,
    ):
        collection = []
        if not zeros_first:
            if n_svd >= 1:
                fitA = self.fit_poles_mod_zeros
            else:
                fitA = self.fit_poles
            if n_svd >= 2:
                fitB = self.fit_zeros_mod_poles
            else:
                fitB = self.fit_zeros
            fitC = self.fit_poles
            fitD = self.fit_zeros
        else:
            if n_svd >= 1:
                fitA = self.fit_zeros_mod_poles
            else:
                fitA = self.fit_zeros
            if n_svd >= 2:
                fitB = self.fit_poles_mod_zeros
            else:
                fitB = self.fit_poles
            fitC = self.fit_zeros
            fitD = self.fit_poles

        fitA(
            max_size  = max_size,
        )
        fitB(
            max_size = max_size,
        )

        for i in range(n_iter - 1):
            fitC()
            fitD()
            if collect_all:
                collection.append(self.copy())

        if n_iter > 0:
            fitC()
            fitD()
        if collect_all:
            collection.append(self.copy())
        return collection

    def fit_pz(
        self,
        n_iter = 0,
        n_svd = 1,
        **kwargs
    ):
        return self.fit_pzpz(
            n_iter = n_iter,
            n_svd  = n_svd,
            **kwargs
        )

    def xfer_eval(self, F_Hz):
        #TODO account for phase difference
        X = np.exp(1j * np.pi * F_Hz / self.F_nyquist_Hz)
        h_val = self.poly.valfromroots(X, self.zeros) / self.poly.valfromroots(X, self.poles)
        return self.gain * h_val


