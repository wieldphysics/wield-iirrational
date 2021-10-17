#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from ...fitters_rational import (
    RationalDiscFilter,
    ChebychevFilter
)

from . import nyquist_move


def rational_disc_fit(
    ZPKrep,
    order           = None,
    order_max       = None,
    order_min       = 20,
    aid             = None,
):
    nyquist_final_Hz = None

    if nyquist_final_Hz is not None and nyquist_final_Hz <= 0:
        raise RuntimeError(
            "nyquist_final_Hz must be positive for Z domain or None"
            " for Sf domain."
        )

    if order is not None:
        order     = int(order)
    order_max = int(order_max)
    order_min = int(order_min)

    N_first = order_min
    N_final = int(min(len(ZPKrep.F_Hz) // 10, order_max))

    if order is not None:
        return ratdisc_single(
            ZPKrep,
            order            = order,
            aid              = aid,
        )

    #otherwise, scan through
    N_current = N_first

    fitter_last = ratdisc_single(
        ZPKrep,
        order           = N_current,
        aid             = aid,
    )
    restot_last = fitter_last.residuals_average

    def fitter_ord(fitter):
        return max(len(fitter.zeros), len(fitter.poles))

    while True:
        if N_current == N_final:
            fitter_use = fitter_last
            break

        N_current = N_current * 2
        if N_current > N_final:
            N_current = N_final

        fitter = ratdisc_single(
            ZPKrep,
            order           = N_current,
            aid             = aid,
        )
        restot = fitter.residuals_average
        fitter_red = fitter.copy()
        fitter_red.matched_pairs_clear(Q_rank_cutoff = .2)
        restot_red = fitter_red.residuals_average

        if restot_last < restot:
            fitter_use = fitter_last
            aid.log("Using last (direct)!", fitter_ord(fitter_last))
            break

        if restot_last < restot_red:
            fitter_use = fitter_last
            aid.log("Using last (reduced)!", fitter_ord(fitter_last))
            break

        if restot_last < 1.10 * restot:
            ord_red = fitter_ord(fitter_red)
            ord_last = fitter_ord(fitter_last)
            if ord_red < ord_last:
                fitter_use = fitter_red
            else:
                fitter_use = fitter_last
            aid.log("Using current")
            break
        #else continue the loop
        fitter_last = fitter
        restot_last = restot

    fitter = fitter_use
    if nyquist_final_Hz == 'keep':
        #this means to not perform the nyquist shift,
        # only really used in testing
        pass
    elif nyquist_final_Hz is None or nyquist_final_Hz == 0:
        ZPKrep = nyquist_move.nyquist_remove(
            fitter,
            split_neg_real = True,
            clear_neg_real = True,
            aid              = aid,
        )
        fitter = ChebychevFilter(
            ZPKrep = ZPKrep,
        )
    elif nyquist_final_Hz > 0:
        fitter = nyquist_move.nyquist_move(
            fitter,
            nyquist_new = nyquist_final_Hz,
            split_neg_real = True,
            clear_neg_real = True,
            aid              = aid,
        )
    else:
        raise RuntimeError((
            "Nyquist cannot be negative. "
            " Positive for Z-domain, None or 0 for S."
            " or \"keep\" to maintain the maxdata nyquist"
        ))

    return fitter


def ratdisc_single(
    ZPKrep,
    order            = 20,
    relative_degree  = 0,
    nyquist_final_Hz = None,
    aid              = None,
):

    if nyquist_final_Hz is not None and nyquist_final_Hz <= 0:
        raise RuntimeError("nyquist_final_Hz must be >0 or None")

    F_nyquist_Hz = 1 * ZPKrep.F_Hz[-1]

    fitter = RationalDiscFilter(
        ZPKrep = ZPKrep,
        npoles        = order,
        nzeros        = order,
        F_nyquist_Hz  = F_nyquist_Hz,
    )
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff = .5)
    fitter.stabilize = True
    #fitter.npoles = len(fitter.poles)
    #fitter.nzeros = len(fitter.zeros)
    ##fitter.stabilize = True
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff = .5)
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    return fitter


