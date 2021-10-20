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
from wavestate import declarative
from .roots_matching import nearest_pairs
from .. import representations


def Q_rank_calc(z, p):
    if p.real == 0 or z.real == 0:
        if p.real == z.real:
            Q_rank = 0
        else:
            #TODO
            #should use the data spacing to regularize this case
            Q_rank = 1e3
    else:
        res_ratio = (z.real/p.real)
        Q_rank = abs(p-z) * (1/(p.real)**2 + 1/(z.real)**2)**.5 + abs(res_ratio - 1/res_ratio)
    return Q_rank


def Q_rank_single(r):
    Q_rank = abs(r.imag / r.real)
    return Q_rank


def order_reduce_zpk(
    zpk,
    Q_rank_cutoff = 1e-5,
    Q_rank_cutoff_unstable = None,
    reduce_c = True,
    reduce_r = False,
    RBalgo = representations.root_bunch.RBalgo,
):
    zeros = RBalgo.expect_atleast(
        zpk[0],
        constraint = RBalgo.root_constraints.mirror_real
    )
    poles = RBalgo.expect_atleast(
        zpk[1],
        constraint = RBalgo.root_constraints.mirror_real
    )
    Pc = poles.c
    Zc = zeros.c
    Pr = poles.r
    Zr = zeros.r

    #print(len(Zc), len(Pc))
    if reduce_c:
        if Q_rank_cutoff_unstable:
            rpB = nearest_pairs(Zc[Zc.real > 0], Pc[Pc.real > 0])
            Zc = list(rpB.l1_remain) + list(Zc[Zc.real <= 0])
            Pc = list(rpB.l2_remain) + list(Pc[Pc.real <= 0])

            removed_rzp_list = []
            for z, p in rpB.r12_list:
                Q_rank = Q_rank_calc(p, z)
                #print(z, p, Q_rank)
                if Q_rank < Q_rank_cutoff_unstable:
                    removed_rzp_list.append((Q_rank, z, p))
                    continue
                Zc.append(z)
                Pc.append(p)
        rpB = nearest_pairs(Zc, Pc)
        Zc = list(rpB.l1_remain)
        Pc = list(rpB.l2_remain)

        removed_rzp_list = []
        for z, p in rpB.r12_list:
            Q_rank = Q_rank_calc(p, z)
            #print(z, p, Q_rank)
            if Q_rank < Q_rank_cutoff:
                removed_rzp_list.append((Q_rank, z, p))
                continue
            Zc.append(z)
            Pc.append(p)

    if reduce_r:
        if Q_rank_cutoff_unstable:
            rpB = nearest_pairs(Zr[Zr.real > 0], Pr[Pr.real > 0])
            Zr = list(rpB.l1_remain) + list(Zr[Zr.real <= 0])
            Pr = list(rpB.l2_remain) + list(Pr[Pr.real <= 0])

            removed_rzp_list = []
            for z, p in rpB.r12_list:
                Q_rank = Q_rank_calc(p, z)
                #print(z, p, Q_rank)
                if Q_rank < Q_rank_cutoff_unstable:
                    removed_rzp_list.append((Q_rank, z, p))
                    continue
                Zr.append(z)
                Pr.append(p)
        rpB = nearest_pairs(Zr, Pr)
        Zr = list(rpB.l1_remain)
        Pr = list(rpB.l2_remain)

        removed_rzp_list = []
        for z, p in rpB.r12_list:
            Q_rank = Q_rank_calc(p, z)
            #print(z, p, Q_rank)
            if Q_rank < Q_rank_cutoff:
                removed_rzp_list.append((Q_rank, z, p))
                continue
            Zr.append(z)
            Pr.append(p)

    return (
        np.concatenate([Zr, Zc, np.conjugate(Zc)]),
        np.concatenate([Pr, Pc, np.conjugate(Pc)]),
        zpk[2]
    )
