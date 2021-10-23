#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Utilities to manipulate ZPK roots S to/from Z, and make transfer functions
"""

import numpy as np
import warnings

from .roots_bin import roots_bin_type

imag_tol = 1e-12


def TF_ZPK(F_Hz, ZPK, F_nyquist_Hz=None, group_len=1):
    """
    Stably generates the ZPK transfer function, interlacing zeros and poles to prevent small number accumulating
    """
    # TODO deprecate
    Z, P, K = ZPK
    if F_nyquist_Hz is None:
        X = 1j * F_Hz
    else:
        X = np.exp(1j * np.pi * F_Hz / F_nyquist_Hz)

    Z = np.asarray(Z)
    P = np.asarray(P)
    mlen = max(len(P), len(Z))
    h = 1
    for idx in range((mlen - 1) // group_len + 1):
        zs = Z[idx * group_len : (idx + 1) * group_len]
        ps = P[idx * group_len : (idx + 1) * group_len]
        h *= np.polynomial.polynomial.polyvalfromroots(X, zs)
        h /= np.polynomial.polynomial.polyvalfromroots(X, ps)
    return K * h


def ZtoS(ZPK, F_nyquist_Hz, F_Hz=None, gain_scale_policy="median", error=False):
    # use a list for mutability due to local/nonlocal referencing semantics
    gain_flip = [1]

    def roots_move(roots_list):
        moved = []
        for r in roots_list:
            if abs(r.imag) < imag_tol:
                if r < 0:
                    if error:
                        raise RuntimeError("Negative real root found in ZtoS")
                    else:
                        warnings.warn("Negative real root found in ZtoS, dropping")
                    gain_flip[0] *= -1
                else:
                    real_Hz = -(1 - r) * F_nyquist_Hz / np.pi
                    moved.append(real_Hz)
            else:
                F_Hz = np.angle(r) / np.pi * F_nyquist_Hz
                real_Hz = -(1 - abs(r)) * F_nyquist_Hz / np.pi
                moved.append(real_Hz + 1j * F_Hz)
        return moved

    Z_z, P_z, K_z = ZPK
    Z_s = roots_move(Z_z)
    P_s = roots_move(P_z)

    if F_Hz is None:
        gain_scale = 1
    else:
        F_Hz = np.asarray(F_Hz)
        xfer_ratio = abs(
            TF_ZPK(F_Hz, (Z_s, P_s, 1), F_nyquist_Hz=None)
            / TF_ZPK(F_Hz, (Z_z, P_z, 1), F_nyquist_Hz=F_nyquist_Hz)
        )
        if gain_scale_policy == "median":
            gain_scale = np.median(xfer_ratio)
        elif gain_scale_policy in ["mean", "average"]:
            gain_scale = np.mean(xfer_ratio)
        else:
            raise RuntimeError(
                'Unrecognized setting for gain_scale_policy, must be "median", "average" (or "Mean").'
            )
    K_s = K_z * gain_flip[0] / gain_scale
    return (np.asarray(Z_s), np.asarray(P_s), K_s)


def StoZ(ZPK, F_nyquist_Hz, F_Hz=None, gain_scale_policy="median", error=False):
    # use a list for mutability due to local/nonlocal referencing semantics
    gain_flip = [1]

    def roots_move(roots_list):
        moved = []
        for r in roots_list:
            if abs(r.imag / r.real) < imag_tol:
                real_D = 1 + r * np.pi / F_nyquist_Hz
                moved.append(real_D)
            else:
                if r.imag >= F_nyquist_Hz:
                    if error:
                        raise RuntimeError(
                            "Root Frequency above the nyquist encountered in StoZ"
                        )
                    else:
                        warnings.warn(
                            "Root Frequency above the nyquist encountered in StoZ, dropping"
                        )
                else:
                    F_D = np.exp(1j * np.pi * r.imag / F_nyquist_Hz)
                    amp_D = 1 + r.real * np.pi / F_nyquist_Hz
                    moved.append(amp_D * F_D)
        return moved

    Z_s, P_s, K_s = ZPK
    Z_z = roots_move(Z_s)
    P_z = roots_move(P_s)
    if F_Hz is None:
        gain_scale = 1
    else:
        F_Hz = np.asarray(F_Hz)
        xfer_ratio = abs(
            TF_ZPK(F_Hz, (Z_z, P_z, 1), F_nyquist_Hz=F_nyquist_Hz)
            / TF_ZPK(F_Hz, (Z_s, P_s, 1), F_nyquist_Hz=None)
        )

        if gain_scale_policy == "median":
            gain_scale = np.median(xfer_ratio)
        elif gain_scale_policy in ["mean", "average"]:
            gain_scale = np.mean(xfer_ratio)
        else:
            raise RuntimeError(
                'Unrecognized setting for gain_scale_policy, must be "median", "average" (or "Mean").'
            )

    K_z = K_s * gain_flip[0] / gain_scale
    return (Z_z, P_z, K_z)


def SorZtoSorZ(
    ZPK,
    F_nyquist_Hz_in,
    F_nyquist_Hz_out,
    F_Hz=None,
    gain_scale_policy="median",
    error=False,
):
    if F_nyquist_Hz_in is not None:
        ZPK_S = ZtoS(
            ZPK,
            F_nyquist_Hz=F_nyquist_Hz_in,
            F_Hz=F_Hz,
            gain_scale_policy=gain_scale_policy,
            error=error,
        )
    else:
        ZPK_S = ZPK

    if F_nyquist_Hz_out is not None:
        ZPK_2 = StoZ(
            ZPK_S,
            F_nyquist_Hz=F_nyquist_Hz_out,
            F_Hz=F_Hz,
            gain_scale_policy=gain_scale_policy,
            error=error,
        )
    else:
        ZPK_2 = ZPK_S
    return ZPK_2


def ZPK_fill(
    ZPK=None,
    Z=None,
    P=None,
    K=None,
    F_nyquist_Hz=None,
):
    if Z is None:
        Z = ZPK[0]
    if P is None:
        P = ZPK[1]
    if K is None:
        K = ZPK[2]

    z_r, z_c, z_u, pol = roots_bin_type(
        Z,
        policy="auto",
        F_nyquist_Hz=F_nyquist_Hz,
    )
    if len(z_u) > 0:
        raise RuntimeError("Found unmatched zeros: {0}".format(str(z_u)))
    zeros = []
    zeros.extend(z_r)
    zeros.extend(z_c)
    zeros.extend([r.conjugate() for r in z_c])

    p_r, p_c, p_u, pol = roots_bin_type(
        P,
        policy="auto",
        F_nyquist_Hz=F_nyquist_Hz,
    )
    if len(p_u) > 0:
        raise RuntimeError("Found unmatched poles: {0}".format(str(p_u)))

    poles = []
    poles.extend(p_r)
    poles.extend(p_c)
    poles.extend([r.conjugate() for r in p_c])
    return [zeros, poles, K]
