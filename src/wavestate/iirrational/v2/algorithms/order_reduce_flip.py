#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from ... import TFmath

from ...utilities import ensure_aid

from . import algorithms


def skip_list_listing(group, non_mindelay=True):
    zr_z = []
    zr_i = []
    for idx_z, z in enumerate(group):
        if z.real < 0 and not non_mindelay:
            continue
        zr_z.append(z)
        zr_i.append(idx_z)
    # print(zr_i, zr_z)
    return zip(zr_i, zr_z)


def ranking_delay_flip(
    aid,
    marginalize_delay=False,
    non_mindelay=True,
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []
    best_worst_list = []

    Zr = aid.fitter.zeros.r
    Zc = aid.fitter.zeros.c

    for idx_z, z in skip_list_listing(Zr, non_mindelay=False):
        z_flip = -z.conjugate()
        pb = algorithms.resrank_program(
            aid.fitter,
            rank_zp_idx_list,
            "flip_r",
            ["ZrDelIdx", idx_z, "ZrAdd", z_flip],
            variant="OrdFl",
            marginalize_delay=marginalize_delay,
            idx=("r", idx_z),
        )
        rank = pb.rank
        best_worst_list.append([rank, "r", idx_z, z, pb])

    for idx_z, z in skip_list_listing(Zc, non_mindelay):
        z_flip = -z.conjugate()
        pb = algorithms.resrank_program(
            aid.fitter,
            rank_zp_idx_list,
            "flip_c",
            ["ZcDelIdx", idx_z, "ZcAdd", z_flip],
            variant="OrdFl",
            marginalize_delay=marginalize_delay,
            idx=("c", idx_z),
        )
        rank = pb.rank
        best_worst_list.append([rank, "c", idx_z, z, pb])

    if not best_worst_list:
        return None

    best_worst_list.sort()
    best_worst_list.reverse()

    N = 0
    while best_worst_list:
        if N > 3:
            break
        bw = best_worst_list.pop()
        rank = bw[0]
        rank_zp_idx_list = bw[4:]
        trials = algorithms.ranking_reduction_trials(
            aid,
            rank_zp_idx_list,
            greedy=True,
            num_total_max=4,
        )
        trial = trials[0]
        did_reduce = aid.fitter_check(
            trial.fitter,
            hint_name="flip",
            variant="OrdC"
        )
        if did_reduce:
            aid.fitter_update(
                trial.fitter,
                representative=True,
            )
            aid.log_progress(
                5,
                ("zero flipped, maxzp {}, residuals={:.2e}, reldeg={}").format(
                    aid.fitter_orders().maxzp, aid.fitter.residuals_average, aid.fitter_orders().reldeg,
                ),
            )
            break
        N += 1

    return did_reduce


def order_reduce_flip(aid, non_mindelay=True):
    aid.log_progress(
        5,
        ("zero flipping, maxzp {}, residuals={:.2e}, {:.2e}, reldeg={}").format(
            aid.fitter_orders().maxzp, aid._fitter_current.residuals_average, aid.fitter.residuals_average, aid.fitter_orders().reldeg,
        ),
    )
    aid.fitter_update(
        representative=True,
    )
    aid.log_progress(
        5,
        ("zero flipping, maxzp {}, residuals={:.2e}, {:.2e}, reldeg={}").format(
            aid.fitter_orders().maxzp, aid._fitter_current.residuals_average, aid.fitter.residuals_average, aid.fitter_orders().reldeg,
        ),
    )
    # import pdb; pdb.set_trace()

    check_sign = aid.fitter.check_sign
    while True:
        ret = ranking_delay_flip(
            aid,
            marginalize_delay=False,
            non_mindelay=non_mindelay,
        )
        if ret is None:
            break
        if not ret:
            break
    aid.fitter.check_sign = check_sign
    return

def order_reduce_flip(aid, non_mindelay=True):
    aid.log_progress(
        5,
        ("zero flipping, maxzp {}, residuals={:.2e}, {:.2e}, reldeg={}").format(
            aid.fitter_orders().maxzp, aid._fitter_current.residuals_average, aid.fitter.residuals_average, aid.fitter_orders().reldeg,
        ),
    )

    from ... import fitters_ZPK
    fitter_new = aid.fitter.regenerate(
        ZPKrep=aid.fitter.ZPKrep,
        coding_map=fitters_ZPK.codings_s.coding_maps.SOSnlpoles,
        check_sign=False,
    )
    with fitter_new.with_codings_only([fitter_new.gain_coding]):
        fitter_new.optimize()
    fitter_new.optimize()
    algorithms.sign_check_flip(fitter_new)
    aid.fitter_update(fitter_new)

    fitter_new = aid.fitter.regenerate(
        ZPKrep=fitter_new.ZPKrep,
        # coding_map=fitters_ZPK.codings_s.coding_maps.SOSnlpoles,
        # check_sign=True,
    )
    aid.fitter_update(
        fitter_new,
        representative=True,
    )
    aid.log_progress(
        5,
        ("zero flipping, maxzp {}, residuals={:.2e}, {:.2e}, reldeg={}").format(
            aid.fitter_orders().maxzp, aid._fitter_current.residuals_average, aid.fitter.residuals_average, aid.fitter_orders().reldeg,
        ),
    )
    return
