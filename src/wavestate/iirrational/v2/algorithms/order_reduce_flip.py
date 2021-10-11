# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

from ... import TFmath

from ..utilities import ensure_aid

from . import algorithms


def skip_list_listing(skip, group):
    zr_z = []
    zr_i = []
    for idx_z, z in enumerate(group):
        if z.real < 0:
            continue
        zr_z.append(z)
        zr_i.append(idx_z)
    for idx in reversed(sorted(TFmath.nearest_idx(skip, zr_z))):
        del zr_z[idx]
        del zr_i[idx]
    #print(zr_i, zr_z)
    return zip(zr_i, zr_z)

def ranking_delay_flip(
    aid,
    marginalize_delay  = True,
    skip_r = [],
    skip_c = [],
):
    aid = ensure_aid(aid)
    rank_zp_idx_list = []
    best_worst_list = []

    Zr = aid.fitter.zeros.r
    Zc = aid.fitter.zeros.c

    for idx_z, z in skip_list_listing(skip_r, Zr):
        z_flip = -z.conjugate()
        pb = algorithms.resrank_program(
            aid.fitter, rank_zp_idx_list,
            'flip_r', ['ZrDel', idx_z, 'ZrAdd', z_flip],
            variant = 'OrdFl',
            marginalize_delay = marginalize_delay,
            idx = ('r', idx_z),
        )
        rank = pb.rank
        best_worst_list.append([rank, 'r', idx_z, z, pb])

    for idx_z, z in skip_list_listing(skip_c, Zc):
        z_flip = -z.conjugate()
        pb = algorithms.resrank_program(
            aid.fitter, rank_zp_idx_list,
            'flip_c', ['ZcDel', idx_z, 'ZcAdd', z_flip],
            variant = 'OrdFl',
            marginalize_delay = marginalize_delay,
            idx = ('c', idx_z),
        )
        rank = pb.rank
        best_worst_list.append([rank, 'c', idx_z, z, pb])

    if not best_worst_list:
        return None

    best_worst_list.sort()
    best_worst_list.reverse()

    aid.log(best_worst_list)
    improved = False
    while best_worst_list:
        bw = best_worst_list.pop()
        rank = bw[0]
        rank_zp_idx_list = bw[4:]
        trials = algorithms.ranking_reduction_trials(
            aid,
            rank_zp_idx_list,
            num_total_max = None,
            num_type_max  = None,
            greedy = True,
        )
        trial = trials[0]
        if trial.improved:
            improved = True
            aid.fitter_update(trial.fitter)
            break
        else:
            if bw[1] == 'c':
                skip_c.append(bw[3])
            else:
                skip_r.append(bw[3])

    return improved, skip_r, skip_c


def order_reduce_flip(aid):
    #TODO, make num_trials a hint value
    skip_r = []
    skip_c = []
    improved = True

    while True:
        ret = ranking_delay_flip(
            aid,
            marginalize_delay  = True,
            skip_r = skip_r,
            skip_c = skip_c,
        )
        if ret is None:
            break
        improved, skip_r, skip_c = ret
