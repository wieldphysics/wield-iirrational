#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from ..utilities import ensure_aid

from .. import TFmath
from .. import representations

from .algorithms import algorithms
from .algorithms import order_reduce_programs


def Q_rank_calc(z, p):
    if p.real == 0 or z.real == 0:
        if p.real == z.real:
            Q_rank = 0
        else:
            # TODO
            # should use the data spacing to regularize this case
            Q_rank = 1e3
    else:
        res_ratio = z.real / p.real
        Q_rank = abs(p - z) * (1 / (p.real) ** 2 + 1 / (z.real) ** 2) ** 0.5 + abs(
            res_ratio - 1 / res_ratio
        )
    return Q_rank


def Q_rank_single(r):
    Q_rank = abs(r.imag / r.real)
    return Q_rank


def order_reduce(
    aid,
    Q_rank_cutoff=1,
    optimize=True,
    reduce_c=True,
    reduce_r=False,
):
    # TODO, can run a high Q_rank_cutoff and the progressively add pairs back in!
    #
    aid = ensure_aid(aid)

    Pc = aid.fitter.poles.c
    Zc = aid.fitter.zeros.c
    Pr = aid.fitter.poles.r
    Zr = aid.fitter.zeros.r

    removed_rzp_list = []
    # print(len(Zc), len(Pc))
    if reduce_c:
        rpB = TFmath.nearest_pairs(Zc, Pc)
        Zc = list(rpB.l1_remain)
        Pc = list(rpB.l2_remain)

        for z, p in rpB.r12_list:
            Q_rank = Q_rank_calc(p, z)
            # print(z, p, Q_rank)
            if Q_rank < Q_rank_cutoff:
                removed_rzp_list.append((Q_rank, z, p))
                continue
            Zc.append(z)
            Pc.append(p)

    if reduce_r:
        rpB = TFmath.nearest_pairs(Zr, Pr)
        Zr = list(rpB.l1_remain)
        Pr = list(rpB.l2_remain)

        for z, p in rpB.r12_list:
            Q_rank = Q_rank_calc(p, z)
            # print(z, p, Q_rank)
            if Q_rank < Q_rank_cutoff:
                removed_rzp_list.append((Q_rank, z, p))
                continue
            Zr.append(z)
            Pr.append(p)

    fitter_new = aid.fitter.regenerate(
        ZPKrep=aid.fitter.ZPKrep,
        z_c=Zc,
        p_c=Pc,
        z_r=Zr,
        p_r=Pr,
    )
    assert(fitter_new.order_relative == aid.fitter.order_relative)
    if optimize:
        with fitter_new.with_codings_only([fitter_new.gain_coding]):
            fitter_new.optimize()
        fitter_new.optimize()
        algorithms.sign_check_flip(fitter_new)
        aid.fitter_update(fitter_new)
    else:
        algorithms.sign_check_flip(fitter_new)
        aid.fitter_update(fitter_new, representative=False)
    
    aid.fitter_checkpoint()
    return removed_rzp_list


def order_restore(
    aid,
    removed_rzp_list,
):
    rank_zp_idx_list = []

    removed_rzp_list.sort()
    for r, z, p in removed_rzp_list:
        algorithms.resrank_program(
            aid.fitter,
            rank_zp_idx_list,
            "zp_restore",
            ["ZcAdd", z, "PcAdd", p],
            marginalize_delay=False,
        )

    did_reduce = True
    ever_restored = False
    while did_reduce:
        trials, rank_zp_idx_list = algorithms.ranking_reduction_trials(
            aid,
            rank_zp_idx_list,
            greedy=True,
            return_remaining=True,
            num_total_max=10,
        )
        if not trials:
            return False
        algorithms.sign_check_flip(trials[0].fitter)
        # print("HMM", aid.fitter.residuals_average, trials[0].fitter.residuals_average)
        trial = trials[0]
        did_reduce = aid.fitter_check(
            trial.fitter,
            hint_name="ordrestore",
            variant=trial.ord_str,
        )
        if did_reduce:
            ever_restored = True
    return ever_restored


def order_reduce_successive(
    aid,
    num_total_max=4,
    num_type_max=2,
    marginalize_delay=True,
    representative=True,
):
    deg_min = aid.hint("total_degree_min")
    if deg_min is None:
        return
    while aid.fitter_orders().maxzp > deg_min:
        greedy_order = aid.hint("greedy_order")
        if greedy_order is None:
            greedy = False
        elif greedy_order < aid.fitter_orders().maxzp:
            greedy = True
        else:
            greedy = False
        rank_zp_idx_list = []

        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_c(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_r(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_c2r(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_cc(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_rr(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_crr(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.sort(key=lambda pb: pb.rank)
        trials = algorithms.ranking_reduction_trials(
            aid,
            rank_zp_idx_list,
            num_total_max=num_total_max,
            num_type_max=num_type_max,
            ranking_factor_max=100,
            greedy=greedy,
        )
        if not trials:
            return False
        trial = trials[0]
        algorithms.sign_check_flip(trial.fitter)
        # print("HMM", aid.fitter.residuals_average, trial.fitter.residuals_average)
        aid.fitter_update(
            trial.fitter,
            representative=representative,
        )
        if not representative:
            aid.fitter_checkpoint()
        aid.log_progress(
            5,
            ("order reduced to {}, residuals={:.2e}").format(
                aid.fitter_orders().maxzp, aid.fitter.residuals_average
            ),
        )
    return


def order_reduce_selective(
    aid,
    num_total_max=4,
    num_type_max=2,
    marginalize_delay=True,
    representative=True,
):
    # TODO, make num_trials a hint value
    ever_reduced = False
    did_reduce = True
    while did_reduce:
        greedy_order = aid.hint("greedy_order")
        if greedy_order is None:
            greedy = False
        elif greedy_order < aid.fitter_orders().maxzp:
            greedy = True
        else:
            greedy = False

        rank_zp_idx_list = []

        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_c(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_r(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_c2r(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_cc(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_rr(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.extend(
            order_reduce_programs.ranking_reduction_crr(
                aid, marginalize_delay=marginalize_delay
            )
        )
        rank_zp_idx_list.sort(key=lambda pb: pb.rank)
        trials = algorithms.ranking_reduction_trials(
            aid,
            rank_zp_idx_list,
            num_total_max=num_total_max,
            num_type_max=num_type_max,
            ranking_factor_max=100,
            greedy=greedy,
        )
        if not trials:
            return False
        trial = trials[0]
        algorithms.sign_check_flip(trial.fitter)
        # print("HMM", aid.fitter.residuals_average, trial.fitter.residuals_average)
        did_reduce = aid.fitter_check(
            trial.fitter,
            hint_name="ordred",
            variant="OrdDn",
        )
        if did_reduce:
            aid.fitter_update(representative=representative)
            aid.fitter_checkpoint()
            aid.log_progress(
                5,
                ("order reduced to {}, residuals={:.2e}, reldeg={}").format(
                    aid.fitter_orders().maxzp, aid.fitter.residuals_average, aid.fitter_orders().reldeg,
                ),
            )
            ever_reduced = True
        else:
            aid.log_progress(6, "order not reduced")
    return ever_reduced
