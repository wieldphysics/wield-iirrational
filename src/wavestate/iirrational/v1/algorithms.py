# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np

from ..utilities import ensure_aid
#from .. import fitters_ZPK


def sign_check_flip(fitter):
    """
    """
    xfer = fitter.xfer_fit
    data = fitter.data
    rat = data / xfer
    rat_ang = np.exp(1j * np.angle(rat))
    ang_avg_fit = np.sum(rat_ang * fitter.W**2) / np.sum(fitter.W**2)

    if (ang_avg_fit.real < 0):
        fitter.gain = -fitter.gain


def flip_mindelay_opt(aid):
    """
    Attempts to flip each non-mindelay zero and then reoptimize
    """
    aid = ensure_aid(aid)

    #TODO, also deal with real roots
    min_idx = 0

    while True:
        coding_lst = list(aid.fitter.num_codings)
        for idx, coding in enumerate(aid.fitter.num_codings):
            if idx <= min_idx:
                #TODO, change logic to not need this and be faster
                continue
            zeros = coding.roots_c_Sf()
            if zeros:
                z = zeros[0]
                if z.real > 0:
                    aid.log("trying flipping", z)
                    #using a coding which cannot flip over,
                    #since that affects the sign of the gain
                    coding_ins = aid.fitter.coding_map.num_c(aid.fitter)
                    #flip the root over, but also reduce its effect to help convergence
                    coding_ins.update_roots_Sf((-z).conjugate())
                    coding_lst[idx] = coding_ins
                    min_idx = idx
                    break
        else:
            #breaks from the while loop only if for-loop break doesn't occur
            break

        fitter_new = aid.fitter.__class__(
            parent = aid.fitter,
            num_codings = coding_lst,
        )
        #TODO, note, should this only be flipped in s-domain?
        #flip the gain,
        #print("GAIN: ", aid.fitter.xfer_fit / fitter_new.xfer_fit)
        #fitter_new.gain = -fitter_new.gain
        #print("GAIN: ", aid.fitter.xfer_fit / fitter_new.xfer_fit)

        with fitter_new.with_codings_only([coding_ins]):
            fitter_new.optimize(aid = aid)
        fitter_new.optimize(aid = aid)
        #print("GAIN3: ", aid.fitter.xfer_fit / fitter_new.xfer_fit)

        aid.fitter_check(
            fitter_new,
            hint_name = 'mindelay_opt',
            variant   = 'OrdC',
        )
    return


def insert_triplets(aid):
    """
    iserts on each poles and zeros a complex and real roots with bandwidth 2x the data
    """
    aid = ensure_aid(aid)

    cplx_t = aid.fitter.coding_map.num_c
    real_t = aid.fitter.coding_map.num_r

    F_l_Hz = aid.fitter.F_max_Hz
    BW_2x_Hz = 2 * F_l_Hz
    F_high_Hz = .90 * F_l_Hz

    coding_num_p = cplx_t(aid.fitter)
    coding_num_p.update_roots_Sf(-BW_2x_Hz + 1j * F_high_Hz)
    coding_num_p2 = real_t(aid.fitter)
    coding_num_p2.update_roots_Sf(-BW_2x_Hz)
    coding_den_p2 = real_t(aid.fitter)
    coding_den_p2.update_roots_Sf(-BW_2x_Hz)
    coding_den_p = cplx_t(aid.fitter)
    coding_den_p.update_roots_Sf(-BW_2x_Hz + 1j * F_high_Hz)

    fitter_new = aid.fitter.__class__(
        parent = aid.fitter,
        num_codings = aid.fitter.num_codings + [coding_num_p2],
        den_codings = aid.fitter.den_codings + [coding_den_p2],
    )
    res_pre = fitter_new.residuals_average
    with fitter_new.with_codings_only([coding_num_p, coding_num_p2] + [coding_den_p, coding_den_p2]):
        fitter_new.optimize(aid = aid)
    res_mid = fitter_new.residuals_average
    fitter_new.optimize(aid = aid)
    res_post = fitter_new.residuals_average

    res_ratio = fitter_new.residuals_average / aid.fitter.residuals_average

    aid.log("TRIPLETS (rat = {0}, pre = {1}, mid = {2}, post = {3})".format(res_ratio, res_pre, res_mid, res_post))
    return aid.fitter_check(
        fitter_new,
        hint_name = 'insert_triplets2',
        variant   = 'OrdUp',
    )


def set_min_BW(aid):
    aid = ensure_aid(aid)

    def modify_codings(codings):
        for coding in codings:
            roots = coding.roots_c_Sf()
            if roots:
                r, = roots
                root_F_Hz = r.imag
                F_idx = np.searchsorted(aid.fitter.F_Hz, root_F_Hz)
                F_idx_low = max(0, F_idx - 2)
                F_idx_high = min(len(aid.fitter.F_Hz) - 1, F_idx + 2)
                #get the local mean width of bins to restrict bandwidths
                F_mean_Hz = (aid.fitter.F_Hz[F_idx_high] - aid.fitter.F_Hz[F_idx_low]) / (F_idx_high - F_idx_low)
                coding.option_set(minimum_BW_Hz = F_mean_Hz * .5)
    modify_codings(aid.fitter.den_codings)
    modify_codings(aid.fitter.num_codings)

    aid.fitter.distance_limit_auto = 2
    aid.fitter_update()
    return
