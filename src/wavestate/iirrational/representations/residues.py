#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Calculation of the RMS of a filter with unit white noise passing through a ZPK.
This is done using an integral and residue calculus. The filter must have more
poles than zeros, unless the ZPK is in the Z domain,
where there is a natural cutoff frequency.
"""

import numpy as np

from scipy.special import factorial
from scipy import signal as scisig
import numpy.polynomial.chebyshev as cheby


class ChebyResiduesHelper(object):
    _ZPK = None

    def __init__(
        self,
        ZPK,
        scale=None,
        rotate=True,
    ):
        # store for the setup, will be cleared after setup has used it
        self._ZPK = ZPK
        self.scale = scale

        if rotate:
            self.rot_p = +1j
            self.rot_n = -1j
        else:
            self.rot_p = 1
            self.rot_n = 1

    def _setup(self):
        if self._ZPK is None:
            return

        zeros, poles, gain = self._ZPK
        # clear the _ZPK holder to indicate that setup has been fully run
        self._ZPK = None
        if self.scale is None:
            self.scale = 1 / max(
                np.amax(np.concatenate([np.abs(zeros), [0]])),
                np.amax(np.concatenate([np.abs(poles), [0]])),
            )
        # print('SCAL:', self.scale)

        self.zeros = self.rot_n * np.asarray(zeros) * self.scale
        self.poles = self.rot_n * np.asarray(poles) * self.scale
        # technically should rescale K here as well, but not, so that it doesn't
        # have to just be inverted later

        b = gain * cheby.chebfromroots(self.zeros)
        self.a = cheby.chebfromroots(self.poles)
        # print("b orig:", cheby.chebfromroots(self.zeros))
        # print("a orig:", self.a)
        self.orderdiff = len(b) - len(self.a)

        k, self.b = cheby.chebdiv(b, self.a)
        self.k = k * 2 ** (self.orderdiff - 1)
        return

    _k_zeros = None

    @property
    def k_zeros(self):
        if self._k_zeros is None:
            self._setup()
            self._k_zeros = self.rot_p * cheby.chebroots(self.k) / self.scale
        return self._k_zeros

    @property
    def k_gain(self):
        self._setup()
        # no need to rescale k here since we didn't rescale it originally
        return self.k[-1]

    def residue_of(self, polemask):
        """
        Takes an untransformed pole, transforms it
        """
        self._setup()
        # rescale pole
        residue_list = []
        bn = self.b.copy()

        an = np.atleast_1d(cheby.chebfromroots(self.poles[polemask]))

        pavg = np.average(self.poles[~polemask])
        # bn(s) / an(s) is (s-po[n])**Nn * b(s) / a(s) where Nn is
        # multiplicity of pole at po[n]
        mult = np.count_nonzero(~polemask)
        for m in range(mult, 0, -1):
            if mult > m:
                # compute next derivative of bn(s) / an(s)
                term1 = cheby.chebmul(cheby.chebder(bn, 1), an)
                term2 = cheby.chebmul(bn, cheby.chebder(an, 1))
                bn = cheby.chebsub(term1, term2)
                an = cheby.chebmul(an, an)
            val = (
                self.rot_p ** (m + self.orderdiff)
                * (
                    (cheby.chebval(pavg, bn) / cheby.chebval(pavg, an))
                    / factorial(mult - m)
                )
                / (self.scale) ** (m + self.orderdiff)
            )

            residue_list.append(val)
        residue_list.reverse()
        return residue_list


def ZPK2residues(
    ZPK,
    F_nyquist_Hz=None,
    Npoles=None,
    _helper_order=None,  # for debugging the derivatives
):
    """
    Npoles allows one to compute residues for fewer poles than given,
    for instance to compute them only for the inner unit circle for ZPKs
    representing power-spectra.
    """
    # currently fills this
    residue_list = []

    zeros, poles, gain = ZPK
    zeros = np.asarray(zeros)
    poles = np.asarray(poles)

    if Npoles is None:
        Npoles = len(poles)

    helper = ChebyResiduesHelper(ZPK)

    # list of enumerated poles that will be depleted as the integral is computed
    # over all residues. Depletion is used since some poles may be stacked
    ppairs = dict(enumerate(poles[:Npoles]))

    res_thresh = 1e-8
    while ppairs:
        idx, p = ppairs.popitem()

        p_eval = p - poles
        pmask = abs(p_eval) > res_thresh
        D_0d = np.prod(p_eval[pmask])
        Nd = np.count_nonzero(~pmask) - 1

        for idx_bad in np.argwhere(~pmask[:Npoles]):
            # remove the multiplicities from the pool to compute
            if idx_bad == idx:
                continue
            # must convert bad to standard python type for the pop to index
            idx_bad = int(idx_bad)
            ppairs.pop(idx_bad)

        if _helper_order is not None and Nd >= _helper_order or Nd >= 4:
            res = helper.residue_of(pmask)
            residue_list.append(
                (p,) + tuple(res),
            )
            # skip the calculation below if already computed
            continue

        N_0d = eval_Nth_derivative(p, zeros, Nd=0)
        R = N_0d / D_0d
        eval_0d = gain * R

        if Nd >= 1:
            # TODO, make the derivatives work up to some order
            N_1d = eval_Nth_derivative(p, zeros, val=N_0d, Nd=1)
            D_1d = eval_Nth_derivative(p, poles, val=D_0d, rootmask=pmask, Nd=1)
            R_1d = (N_1d - D_1d * N_0d / D_0d) / D_0d
            eval_1d = gain * R_1d

        if Nd >= 2:
            N_2d = eval_Nth_derivative(p, zeros, val=N_0d, Nd=2)
            D_2d = eval_Nth_derivative(p, poles, val=D_0d, rootmask=pmask, Nd=2)
            # TODO, make the derivatives work up to some order
            R_2d = (N_2d - D_2d * N_0d / D_0d - 2 * R_1d * D_1d) / D_0d
            eval_2d = gain * R_2d / 2

        if Nd >= 3:
            # TODO, make the derivatives work up to some order
            N_3d = eval_Nth_derivative(p, zeros, val=N_0d, Nd=3)
            D_3d = eval_Nth_derivative(p, poles, val=D_0d, rootmask=pmask, Nd=3)
            # TODO, make the derivatives work up to some order
            R_3d = (
                N_3d
                + (
                    -D_1d * N_2d
                    - D_2d * N_1d
                    - D_3d * N_1d
                    + 2 * R_1d * D_1d ** 2
                    + 2 * N_0d * D_2d * D_1d / D_0d
                )
                / D_0d
                - 2 * (R_2d * D_1d + R_1d * D_2d)
            ) / D_0d
            eval_3d = gain * R_3d / 6

        if Nd == 0:
            residue_list.append(
                (p, eval_0d),
            )
        elif Nd == 1:
            residue_list.append(
                (p, eval_1d, eval_0d),
            )
        elif Nd == 2:
            residue_list.append(
                (p, eval_2d, eval_1d, eval_0d),
            )
        elif Nd == 3:
            residue_list.append(
                (p, eval_3d, eval_2d, eval_1d, eval_0d),
            )
        else:
            raise RuntimeError(
                "Residues with poles of mult>3 not implemented"
                " with direct derivatives"
            )

    if len(zeros) < len(poles):
        return residue_list, ((), (), 0)
    elif len(zeros) == len(poles):
        return residue_list, ((), (), gain)
    else:
        return residue_list, (
            tuple(helper.k_zeros),
            (),
            np.asarray(helper.k_gain).item(),
        )


def eval_Nth_derivative(z, roots, rootmask=None, val=None, Nd=0):
    if rootmask is None:
        rootmask = np.ones(roots.shape, dtype=bool)

    if Nd == 0:
        if val is None:
            return np.prod(z - roots[rootmask])
        else:
            return val

    if val is None and Nd < 3:
        r_eval = z - roots[rootmask]
        val = np.prod(r_eval)

    # TODO, currently this algorithm is very inefficient
    # for large derivatives, it should be able to factor in
    # multiplicities by only iterating indices ABOVE the first zero mask.
    d_sum = 0
    for idx in np.argwhere(rootmask):
        mask_down = np.copy(rootmask)
        mask_down[idx] = 0
        if val is not None:
            r_eval = z - roots[idx]
            # check numerical stability, don't divide out small numbers
            if abs(r_eval) > 1e-8:
                val_down = val / r_eval
            else:
                val_down = None
        else:
            val_down = None
        d_sum += eval_Nth_derivative(z, roots, mask_down, val=val_down, Nd=Nd - 1)
    return np.asarray(d_sum).item()


def ZPK2residues_scipy(
    ZPK,
    F_nyquist_Hz=None,
):
    """
    Useful for testing
    """
    residue_list = []
    b, a = scisig.zpk2tf(*ZPK)
    if F_nyquist_Hz is None:
        r, p, k = scisig.residue(b, a)
    else:
        r, p, k = scisig.residuez(b, a)

    # some indexing tricks to change the output to the form above
    same = p[1:] == p[:-1]
    same_adj = np.concatenate([[False], same, [False]])
    same_shift = ~(same_adj[1:] | same_adj[:-1])
    same_starts = ~same_adj[:-1] & same_adj[1:]
    same_ends = same_adj[:-1] & ~same_adj[1:]
    # the addition of same_adj puts all of the non-multiple back in
    for idx_start, idx_end in zip(
        np.argwhere(same_starts | same_shift), np.argwhere(same_ends | same_shift)
    ):
        idx_start = int(idx_start)
        idx_end = int(idx_end)
        pole = p[idx_start]
        res = r[idx_start : idx_end + 1]
        residue_list.append((pole,) + tuple(res))

    if len(k) < 2:
        return residue_list, ((), (), k[0])
    else:
        Z, P, K = scisig.tf2zpk(k, [1])
        return residue_list, (tuple(Z), tuple(P), K)


def residues2ZPK(rlist, ZPKadd=None, scale=None, rotate=True):
    """ """
    poles = []
    for rtup in rlist:
        pole = rtup[0]
        poles.extend([pole] * (len(rtup) - 1))
    poles = np.asarray(poles)

    if scale is None:
        scale = 1 / np.amax(np.concatenate([abs(poles), [0]]))

    if rotate:
        rot_p = +1j
        rot_n = -1j
    else:
        rot_p = 1
        rot_n = 1

    # scale = .5
    # print(scale)
    poles_scaled = rot_n * poles * scale

    if ZPKadd is None:
        b = [0]
        orderdiff = 1 - len(poles)
    else:
        Z, P, K = ZPKadd
        Z = rot_n * np.asarray(Z) * scale
        b = cheby.chebmul(cheby.chebfromroots(poles_scaled), cheby.chebfromroots(Z) * K)

        orderdiff = len(b) - len(poles) - 2
        if len(P):
            raise NotImplementedError("Currently Cannot Sum ZPK with poles")
    idx_plist = 0
    # print("orderdiff", orderdiff)
    # print("b: start ", b)
    for rtup in rlist:
        pole = poles_scaled[idx_plist]
        mult = len(rtup) - 1
        poles_temp = np.concatenate(
            [
                poles_scaled[:idx_plist],
                poles_scaled[idx_plist + mult :],
            ]
        )
        a_outer = cheby.chebfromroots(poles_temp)

        for m, res in enumerate(rtup[1:]):
            a_inner = cheby.chebfromroots([pole] * (mult - m - 1))
            a_temp = cheby.chebmul(a_outer, a_inner)
            # TODO, might be able to avoid the scaling factor here if scaling
            # of the original b is done differently at the beginning
            res_scaled = res * rot_n ** (m - orderdiff) * scale ** (2 + orderdiff + m)
            b = cheby.chebadd(b, res_scaled * a_temp)
            # print(res_scaled * a_temp)
        idx_plist += mult

    # print(b)
    gain = b[-1] * 2 ** (len(b) - 2)  # / scale**len(b)
    zeros = rot_p * cheby.chebroots(b) / scale

    return (tuple(zeros), tuple(poles), gain)
