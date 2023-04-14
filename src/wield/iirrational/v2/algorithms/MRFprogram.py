#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
import numpy as np
from wield.bunch import Bunch

from ... import TFmath
from ... import representations

from ...utilities import ensure_aid
from ...fitters_ZPK.ZPKrep2MRF import ZP2codings


def resrank_program(
    fitter,
    rank_zp_idx_list,
    name,
    program,
    variant=None,
    marginalize_delay=True,
    **kwargs
):

    Pc = fitter.poles.c
    Zc = fitter.zeros.c
    Pr = fitter.poles.r
    Zr = fitter.zeros.r

    Pc_mod = []
    Pr_mod = []
    Zc_mod = []
    Zr_mod = []

    program_copy = list(program)
    while program_copy:
        arg = program_copy.pop()
        which = program_copy.pop()
        if which == "PcDelIdx":
            Pc_mod.append(Pc[arg])
        elif which == "ZcDelIdx":
            Zc_mod.append(Zc[arg])
        elif which == "PrDelIdx":
            Pr_mod.append(Pr[arg])
        elif which == "ZrDelIdx":
            Zr_mod.append(Zr[arg])
        # swap on these since they are added
        elif which == "PrAdd":
            Zr_mod.append(arg)
        elif which == "ZrAdd":
            Pr_mod.append(arg)
        elif which == "PcAdd":
            Zc_mod.append(arg)
        elif which == "ZcAdd":
            Pc_mod.append(arg)

    PrB = representations.asMRRB(c=Pc_mod, r=Pr_mod)
    ZrB = representations.asMRRB(c=Zc_mod, r=Zr_mod)
    h, lnG = PrB.val_lnG(fitter.F_Hz)
    h, lnG = ZrB.val_lnG(fitter.F_Hz, h=1 / h, lnG=-lnG)

    R = fitter.xfer_fit / (fitter.data * h)

    # this is the exact solution for the gain adjustment required, using the DualB
    # residuals
    G = (
        np.sum(fitter.W ** 2 / TFmath.abs_sq(R))
        / np.sum(fitter.W ** 2 * TFmath.abs_sq(R))
    ) ** 0.25
    R = G * R

    rank = np.sum(TFmath.abs_sq(fitter.residuals_NLmap(R, W=fitter.W)))
    if marginalize_delay:
        # linear detrends the phasing term to remove delay effects
        f = fitter.F_Hz
        # the factor of 2 is not with the exact solution, but seems to work better
        R.imag -= (
            f * np.sum(f * R.imag * fitter.W ** 2) / np.sum(f ** 2 * fitter.W ** 2) / 2
        )
        rank_moddelay = np.sum(TFmath.abs_sq(fitter.residuals_NLmap(R, W=fitter.W)))
        rank = min(rank, rank_moddelay)
    pbunch = Bunch(**kwargs)
    pbunch.rank = rank
    pbunch.name = name
    pbunch.variant = variant
    pbunch.program = program
    pbunch.flip = True
    rank_zp_idx_list.append(pbunch)
    return pbunch


def ranking_reduction_trials(
    aid,
    rank_zp_idx_list,
    num_total_max=4,
    num_type_max=2,
    num_try2=2,
    ranking_factor_max=None,
    greedy=False,
    return_remaining=False,
    reset_delay=True,
    hint_name=None,
):
    aid = ensure_aid(aid)
    if not rank_zp_idx_list:
        if return_remaining:
            return [], rank_zp_idx_list
        else:
            return []

    type_count = dict()
    num_total = 0
    idx_total = 0

    trials = []
    trials_try1 = []
    trials_try2 = []
    ranks_original = []
    minranking = float("inf")

    aid.log_progress(6, "trials started")
    for pbunch in rank_zp_idx_list:
        idx_total += 1

        if ranking_factor_max is not None and (
            pbunch.rank > minranking * ranking_factor_max
        ):
            continue
        if pbunch.rank < minranking:
            minranking = pbunch.rank

        count = type_count.get(pbunch.name, 0) + 1
        if num_type_max is not None and count > num_type_max:
            continue
        type_count[pbunch.name] = count
        num_total += 1
        if num_total_max is not None and num_total > num_total_max:
            break
        ranks_original.append(pbunch.rank)

        Pc_mod = list(aid.fitter.poles.c)
        Zc_mod = list(aid.fitter.zeros.c)
        Pr_mod = list(aid.fitter.poles.r)
        Zr_mod = list(aid.fitter.zeros.r)
        Pc_idx_remove = set()
        Zc_idx_remove = set()
        Pr_idx_remove = set()
        Zr_idx_remove = set()
        Pc_new = []
        Zc_new = []
        Pr_new = []
        Zr_new = []
        # work backward on the command sequence removing the idxs
        trial = Bunch()
        trial.prog_redux = []
        trial.ord_ch = 0
        trial.pbunch = pbunch
        program = list(pbunch.program)
        while program:
            arg = program.pop()
            which = program.pop()
            trial.prog_redux.append(which)
            if which == "PcDelIdx":
                trial.prog_redux.append(Pc_mod[arg])
                trial.ord_ch -= 2
                Pc_idx_remove.add(arg)
            elif which == "ZcDelIdx":
                trial.prog_redux.append(Zc_mod[arg])
                trial.ord_ch -= 2
                Zc_idx_remove.add(arg)
            elif which == "PrDelIdx":
                trial.prog_redux.append(Pr_mod[arg])
                trial.ord_ch -= 1
                Pr_idx_remove.add(arg)
            elif which == "ZrDelIdx":
                trial.prog_redux.append(Zr_mod[arg])
                trial.ord_ch -= 1
                Zr_idx_remove.add(arg)
            elif which == "PrAdd":
                trial.prog_redux.append(arg)
                trial.ord_ch += 1
                Pr_new.append(arg)
            elif which == "ZrAdd":
                trial.prog_redux.append(arg)
                trial.ord_ch += 1
                Zr_new.append(arg)
            elif which == "PcAdd":
                trial.prog_redux.append(arg)
                trial.ord_ch += 2
                Pc_new.append(arg)
            elif which == "ZcAdd":
                trial.prog_redux.append(arg)
                trial.ord_ch += 2
                Zc_new.append(arg)
        for idx in reversed(sorted(Pc_idx_remove)):
            del Pc_mod[idx]
        for idx in reversed(sorted(Zc_idx_remove)):
            del Zc_mod[idx]
        for idx in reversed(sorted(Pr_idx_remove)):
            del Pr_mod[idx]
        for idx in reversed(sorted(Zr_idx_remove)):
            del Zr_mod[idx]

        if trial.ord_ch == 0:
            trial.ord_str = "OrdC"
        elif trial.ord_ch < 0:
            trial.ord_str = "OrdDn"
        else:
            trial.ord_str = "OrdUp"

        coding_map, num_codings_mod, den_codings_mod = ZP2codings(
            aid.fitter,
            zeros=representations.asMRRB(r=Zr_mod, c=Zc_mod),
            poles=representations.asMRRB(r=Pr_mod, c=Pc_mod),
        )
        coding_map, num_codings_new, den_codings_new = ZP2codings(
            aid.fitter,
            zeros=representations.asMRRB(r=Zr_new, c=Zc_new),
            poles=representations.asMRRB(r=Pr_new, c=Pc_new),
            coding_map=coding_map,
        )
        fitter = coding_map.mrf_default(
            parent=aid.fitter,
            num_codings=num_codings_mod + num_codings_new,
            den_codings=den_codings_mod + den_codings_new,
        )

        trial.fitter = fitter
        trial.codings_mod = num_codings_mod + den_codings_new
        trial.codings_new = num_codings_mod + den_codings_new
        trial.rank_original = pbunch.rank
        trials_try1.append(trial)

    def trial_optimize(trial):
        improved = False
        fitter = trial.fitter
        try:
            with fitter.with_codings_only([fitter.gain_coding]):
                fitter.optimize()
        except Exception as e:
            aid.log_debug(9, "Optimize Exception", e)

        try:
            # anneal by moving the original codings first
            #if trial.codings_new:
            #    with fitter.with_codings_only(trial.codings_mod):
            #        fitter.optimize()
            if reset_delay:
                aid.fitter.delay_s = aid.hint("delay_s")
            #fitter.optimize_NM()
            fitter.optimize()
        except Exception as e:
            aid.log_debug(9, "Optimize Exception", e)
            trial.improved = False
            return trial

        if trial.pbunch.variant is not None:
            ord_str = trial.pbunch.variant
        else:
            ord_str = trial.ord_str

        improved = aid.fitter_check(
            fitter,
            variant=ord_str,
            update=False,
            validate=False,
            hint_name=hint_name,
        )
        trial.improved = improved
        return trial

    mt = aid.hint("multithreading", None)
    if mt is not None and mt > 1:
        import multiprocessing.pool

        pool = multiprocessing.pool.ThreadPool(
            processes=mt,
        )
        trial_map = pool.imap_unordered
    else:
        trial_map = map

    improved = False
    for trial in trial_map(trial_optimize, trials_try1):
        aid.log(8, "Reducing Program", pbunch.rank, trial.prog_redux)
        if not trial.improved:
            trials_try2.append(trial)
        else:
            improved = True
            trials_try2.append(trial)
        trials.append(trial)
        if greedy and trial.improved:
            break

    aid.log_progress(7, "trials annealing")

    def trial_optimize_anneal(trial):
        fitter = trial.fitter
        # TODO fix
        from . import algorithms

        if not algorithms.optimize_anneal(aid, fitter):
            return trial
        improved = aid.fitter_check(
            fitter,
            variant=trial.ord_str,
            update=False,
            validate=False,
        )
        trial.improved = improved
        return trial

    if not (greedy and improved):
        trials_try2.sort(key=lambda trial: trial.fitter.residuals_average)
        for trial in trial_map(trial_optimize_anneal, trials_try2[:num_try2]):
            if greedy and trial.improved:
                break

    trials.sort(key=lambda trial: trial.fitter.residuals_average)
    if not return_remaining:
        return trials
    else:
        return trials, rank_zp_idx_list[idx_total:]
