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

from ... import TFmath
from ...utilities import ensure_aid

from . import algorithms


def Q_rank_calc(z, p):
    if p.real == 0 or z.real == 0:
        if p.real == z.real:
            Q_rank = 0
        else:
            #TODO
            #should use the data spacing to regularize this case
            Q_rank = 1e3
    else:
        Q_rank = abs(p-z) * (1/(p.real)**2 + 1/(z.real)**2)**.5
    return Q_rank


def Q_rank_single(r):
    Q_rank = abs(r.imag / r.real)
    return Q_rank


def ranking_reduction_cc(
    aid,
    marginalize_delay = True,
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []

    Pc = aid.fitter.poles.c
    Zc = aid.fitter.zeros.c

    #can always works since this is a balanced reduction

    #print(len(Zc), len(Pc))
    pairs = set()
    lZP = TFmath.nearest_idx(Zc, Pc)
    pairs.update(enumerate(lZP))

    lPZ = TFmath.nearest_idx(Pc, Zc)
    for idx_p, idx_z in enumerate(lPZ):
        pairs.add((idx_z, idx_p))

    lPZ = TFmath.nearest_idx(Pc, -Zc.conjugate())
    #only include ones which are pairing unstable with stable
    for idx_p, idx_z in enumerate(lPZ):
        if Pc[idx_p].real * Zc[idx_z].real < 0:
            pairs.add((idx_z, idx_p))

    for idx_z, idx_p in pairs:
        if idx_z is None or idx_p is None:
            continue
        algorithms.resrank_program(
            aid.fitter, rank_zp_idx_list,
            'cc', ['ZcDelIdx', idx_z, 'PcDelIdx', idx_p],
            marginalize_delay = marginalize_delay
        )
    return rank_zp_idx_list


def ranking_reduction_c(
    aid,
    marginalize_delay = True,
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []

    Pc = aid.fitter.poles.c
    Zc = aid.fitter.zeros.c

    orders = aid.fitter_orders()
    rdmax = aid.hint('relative_degree_max')
    rdmin = aid.hint('relative_degree_min')

    if rdmin is None or orders.reldeg - 2 >= rdmin:
        for idx_z, z in enumerate(Zc):
            algorithms.resrank_program(
                aid.fitter, rank_zp_idx_list,
                'c', ['ZcDelIdx', idx_z],
                marginalize_delay = marginalize_delay
            )

    if rdmax is None or orders.reldeg + 2 <= rdmax:
        for idx_p, p in enumerate(Pc):
            algorithms.resrank_program(
                aid.fitter, rank_zp_idx_list,
                'c', ['PcDelIdx', idx_p],
                marginalize_delay = marginalize_delay
            )
    return rank_zp_idx_list


def ranking_reduction_r(
    aid,
    marginalize_delay = True,
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []

    Pr = aid.fitter.poles.r
    Zr = aid.fitter.zeros.r

    orders = aid.fitter_orders()
    rdmax = aid.hint('relative_degree_max')
    rdmin = aid.hint('relative_degree_min')

    if rdmin is None or orders.reldeg - 1 >= rdmin:
        for idx_z, z in enumerate(Zr):
            algorithms.resrank_program(
                aid.fitter, rank_zp_idx_list,
                'r', ['ZrDelIdx', idx_z],
                marginalize_delay = marginalize_delay
            )

    if rdmax is None or orders.reldeg + 1 <= rdmax:
        for idx_p, p in enumerate(Pr):
            algorithms.resrank_program(
                aid.fitter, rank_zp_idx_list,
                'r', ['PrDelIdx', idx_p],
                marginalize_delay = marginalize_delay
            )
    return rank_zp_idx_list


def ranking_reduction_rr(
    aid,
    marginalize_delay = True,
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []

    Pr = aid.fitter.poles.r
    Zr = aid.fitter.zeros.r

    #doesn't need to check reldeg since this preserves it

    pairs = set()
    lZP = TFmath.nearest_idx(Zr, Pr)
    pairs.update(enumerate(lZP))

    lPZ = TFmath.nearest_idx(Pr, Zr)
    for idx_p, idx_z in enumerate(lPZ):
        pairs.add((idx_z, idx_p))

    #TODO, include these via an option in case of starvation?
    lPZ = TFmath.nearest_idx(Pr, -Zr)
    #only include ones which are pairing unstable with stable
    for idx_p, idx_z in enumerate(lPZ):
        if Pr[idx_p] * Zr[idx_z] < 0:
            pairs.add((idx_z, idx_p))

    for idx_z, idx_p in pairs:
        if idx_z is None or idx_p is None:
            continue
        algorithms.resrank_program(
            aid.fitter, rank_zp_idx_list,
            'rr', ['ZrDelIdx', idx_z, 'PrDelIdx', idx_p],
            marginalize_delay = marginalize_delay
        )
    return rank_zp_idx_list


def ranking_reduction_c2r(
    aid,
    marginalize_delay  = True,
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []

    Pc = aid.fitter.poles.c
    Zc = aid.fitter.zeros.c

    orders = aid.fitter_orders()
    rdmax = aid.hint('relative_degree_max')
    rdmin = aid.hint('relative_degree_min')

    if rdmin is None or orders.reldeg - 1 >= rdmin:
        for idx_z, z in enumerate(Zc):
            Q_rank = Q_rank_single(z)
            if Q_rank > 2:
                continue

            z_r = -abs(z)
            algorithms.resrank_program(
                aid.fitter, rank_zp_idx_list,
                'c2r', ['ZcDelIdx', idx_z, 'ZrAdd', z_r],
                marginalize_delay = marginalize_delay
            )

    if rdmax is None or orders.reldeg + 1 <= rdmax:
        for idx_p, p in enumerate(Pc):
            Q_rank = Q_rank_single(p)
            if Q_rank > 2:
                continue

            p_r = -abs(p)
            algorithms.resrank_program(
                aid.fitter, rank_zp_idx_list,
                'c2r', ['PcDelIdx', idx_p, 'PrAdd', p_r],
                marginalize_delay = marginalize_delay
            )

    return rank_zp_idx_list


def ranking_reduction_crr(
    aid,
    marginalize_delay  = True,
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []

    Pc = aid.fitter.poles.c
    Zc = aid.fitter.zeros.c
    Pr = aid.fitter.poles.r
    Zr = aid.fitter.zeros.r

    Zc_argQsort = np.argsort(Q_rank_single(Zc))
    Pc_argQsort = np.argsort(Q_rank_single(Pc))

    Zr_BWsort = abs(Zr.real)
    Zr_argBWsort = np.argsort(Zr_BWsort)
    Zr_BWsort = Zr_BWsort[Zr_argBWsort]

    Pr_BWsort = abs(Pr.real)
    Pr_argBWsort = np.argsort(Pr_BWsort)
    Pr_BWsort = Pr_BWsort[Pr_argBWsort]

    def crr_ranker(
        Xc, Xc_argQsort,
        Yr, Yr_BWsort, Yr_argBWsort,
        Xc_code, Yr_code, include_unbalanced,
    ):
        if len(Yr_BWsort) <= 1:
            return
        for idx_x in Xc_argQsort[:3]:
            x = Xc[idx_x]
            bw_eff = abs(x)
            y_sort_idx = np.searchsorted(Yr_BWsort, bw_eff)

            if len(Yr_BWsort) >= 2:
                if y_sort_idx == 0:
                    idx_y1 = 0
                    idx_y2 = 1
                elif y_sort_idx >= len(Yr_BWsort) - 1:
                    idx_y1 = len(Yr_BWsort) - 1
                    idx_y2 = len(Yr_BWsort) - 2
                else:
                    idx_y1 = y_sort_idx
                    idx_y2 = y_sort_idx + 1
                #remap from sorted idx back into the original list idx
                idx_y1 = Yr_argBWsort[idx_y1]
                idx_y2 = Yr_argBWsort[idx_y2]
                #ensure that the ordering of the program is guaranteed
                if idx_y1 > idx_y2:
                    idx_y1, idx_y2 = idx_y2, idx_y1

                algorithms.resrank_program(
                    aid.fitter, rank_zp_idx_list,
                    'crr', [Xc_code, idx_x, Yr_code, idx_y1, Yr_code, idx_y2],
                    marginalize_delay = marginalize_delay
                )

            if include_unbalanced and len(Yr_BWsort) >= 1:
                if y_sort_idx == 0:
                    idx_y = 0
                elif y_sort_idx >= len(Yr_BWsort) - 1:
                    idx_y = len(Yr_BWsort) - 1
                else:
                    if bw_eff - Yr[y_sort_idx] < Yr[y_sort_idx + 1] - bw_eff:
                        idx_y = y_sort_idx
                    else:
                        idx_y = y_sort_idx + 1
                #remap from sorted idx back into the original list idx
                idx_y = Yr_argBWsort[idx_y]

                algorithms.resrank_program(
                    aid.fitter, rank_zp_idx_list,
                    'cr', [Xc_code, idx_x, Yr_code, idx_y],
                    marginalize_delay = marginalize_delay
                )

    orders = aid.fitter_orders()
    rdmax = aid.hint('relative_degree_max')
    rdmin = aid.hint('relative_degree_min')

    crr_ranker(
        Xc = Zc, Xc_argQsort = Zc_argQsort,
        Yr = Pr, Yr_BWsort = Pr_BWsort, Yr_argBWsort = Pr_argBWsort,
        Xc_code = 'ZcDelIdx', Yr_code = 'PrDelIdx',
        include_unbalanced = rdmin is None or orders.reldeg - 1 >= rdmin,
    )

    crr_ranker(
        Xc = Pc, Xc_argQsort = Pc_argQsort,
        Yr = Zr, Yr_BWsort = Zr_BWsort, Yr_argBWsort = Zr_argBWsort,
        Xc_code = 'PcDelIdx', Yr_code = 'ZrDelIdx',
        include_unbalanced = rdmax is None or orders.reldeg + 1 <= rdmax,
    )

    return rank_zp_idx_list

