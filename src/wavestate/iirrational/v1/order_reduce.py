# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
import copy
import declarative

from ..utilities import ensure_aid

from .. import fitters_ZPK
from .. import TFmath
from .. import plots
from .. import annotate
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

def match_pairs(
    aid,
    duals = True,
):
    """
    Match unique closest pairs, if they are within a bandwidth of 0Hz, then they are ignored
    """
    aid = ensure_aid(aid)

    def coding_to_roots(coding_lst):
        coding_roots = []
        for coding in coding_lst:
            r_lst = coding.roots_c_Sf()
            if len(r_lst) == 0:
                coding_roots.append(None)
                continue
            elif len(r_lst) == 1:
                r = r_lst[0]
                if r.real > 0:
                    #TODO allow for unstable
                    #coding_roots.append(1/r[0].conjugate())
                    coding_roots.append(None)
                    continue
                else:
                    #then must be within 45deg
                    if r.imag < r.real:
                        #then this root is within a bandwidth of zero and is effectively a real root
                        #these require much more careful consideration
                        #TODO, check that the usage of roots_Sf is OK
                        coding_roots.append(None)
                        continue
                    coding_roots.append(r)
                    continue
            else:
                #TODO account for
                coding_roots.append(None)
                assert(False)
                continue
        return coding_roots

    coding_poles = coding_to_roots(aid.fitter.den_codings)
    coding_zeros = coding_to_roots(aid.fitter.num_codings)

    def nearest_idx(lst_1, lst_2):
        nearest_lst = []
        for r1 in lst_1:
            if r1 is None:
                nearest_lst.append(None)
                continue
            dist_nearest = float('inf')
            idx_nearest = None
            for idx_2, r2 in enumerate(lst_2):
                if r2 is None:
                    continue
                dist = abs(r1 - r2)
                if dist < dist_nearest:
                    idx_nearest = idx_2
                    dist_nearest = dist
            nearest_lst.append(idx_nearest)
        return nearest_lst
    z_nearest = nearest_idx(coding_zeros, coding_poles)
    p_nearest = nearest_idx(coding_poles, coding_zeros)

    if duals:
        #TODO, refactor the duals generation
        duals = []
        for idx_z, idx_p in enumerate(z_nearest):
            if idx_p is None:
                continue
            #coding_z = aid.fitter.num_codings[idx_z]
            #coding_p = aid.fitter.den_codings[idx_p]
            #TODO annotate about stability
            p = coding_poles[idx_p]
            z = coding_zeros[idx_z]
            Q_rank = Q_rank_calc(z, p)
            Q_rank_BW = (TFmath.abs_sq(p-z) + (p.real)**2 + (p.imag)**2)**.5
            F_loc_Hz = (p.imag + z.imag) / 2
            if idx_z == p_nearest[idx_p]:
                duals.append(declarative.Bunch(
                    idx_z        = idx_z,
                    idx_p        = idx_p,
                    Q_rank       = Q_rank,
                    Q_rank_BW    = Q_rank_BW,
                    F_loc_Hz     = F_loc_Hz,
                    dist         = abs(p-z),
                    z            = z,
                    p            = p,
                ))
        duals.sort(key = lambda b : b.Q_rank)
        return duals
    else:
        pairs = []
        def loop(idx_z, idx_p):
            p = coding_poles[idx_p]
            z = coding_zeros[idx_z]
            Q_rank = Q_rank_calc(z, p)
            Q_rank_BW = (TFmath.abs_sq(p-z) + (p.real)**2 + (p.imag)**2)**.5
            F_loc_Hz = (p.imag + z.imag) / 2
            pairs.append(declarative.Bunch(
                idx_z        = idx_z,
                idx_p        = idx_p,
                Q_rank       = Q_rank,
                Q_rank_BW    = Q_rank_BW,
                F_loc_Hz     = F_loc_Hz,
                dist         = abs(p-z),
                z            = z,
                p            = p,
            ))
        for idx_z, idx_p in enumerate(z_nearest):
            if idx_p is None:
                continue
            loop(idx_z, idx_p)
        for idx_p, idx_z in enumerate(p_nearest):
            if idx_z is None:
                continue
            if idx_p == z_nearest[idx_z]:
                #would be redundant
                continue
            loop(idx_z, idx_p)
        pairs.sort(key = lambda b : b.Q_rank)
        return pairs


def add_resrank(
    dual_b,
    aid,
    relupdown = 2,
):
    aid = ensure_aid(aid)
    dual_b.idx_min = np.searchsorted(aid.fitter.F_Hz, dual_b.F_loc_Hz - relupdown * dual_b.Q_rank_BW)
    dual_b.idx_max = np.searchsorted(aid.fitter.F_Hz, dual_b.F_loc_Hz + relupdown * dual_b.Q_rank_BW)

    #aid.log(idx_min, idx_max)
    select = slice(dual_b.idx_min, dual_b.idx_max)
    #F_Hz = aid.fitter.F_Hz[select]
    Xsf   = 1j * aid.fitter.F_Hz[select]
    W    = aid.fitter.W[select]
    res  = aid.fitter.residuals_preferred[select]
    h    = (
        TFmath.polyvalfromroots(Xsf, [dual_b.z])
        / TFmath.polyvalfromroots(Xsf, [dual_b.p])
    )
    #aid.log(abs(h))
    res_h = W * (h.conjugate()**-2 - h)
    R = aid.fitter.xfer_fit[select]/aid.fitter.data[select]
    res = W * ((R.conjugate())**-1 - R)
    R_mod = aid.fitter.xfer_fit[select]/h/aid.fitter.data[select]
    res_mod = W * (R_mod.conjugate()**-1 - R_mod)

    dual_b.res_tot = np.sum(TFmath.abs_sq(res)) / (dual_b.idx_max - dual_b.idx_min)
    dual_b.res_h_tot = np.sum(TFmath.abs_sq(res_h)) / (dual_b.idx_max - dual_b.idx_min)
    dual_b.res_mod_tot = np.sum(TFmath.abs_sq(res_mod)) / (dual_b.idx_max - dual_b.idx_min)
    return dual_b


def order_reduce(
    aid,
    duals = None,
    Q_rank_cutoff = 1,
    rolloff_order = 0,
):
    aid = ensure_aid(aid)

    if duals is None:
        duals = match_pairs(aid.fitter)
    r1 = TFmath.abs_sq(aid.fitter.residuals.resP)

    num_codings       = []
    den_codings       = []
    exclude_p = set()
    exclude_z = set()
    duals_used = []
    duals_rejected = []
    for b in duals:
        used = False
        if b.Q_rank < Q_rank_cutoff:
            exclude_p.add(b.idx_p)
            exclude_z.add(b.idx_z)
            duals_used.append(b)
            used = True
        else:
            F_idx = np.searchsorted(aid.fitter.F_Hz, b.F_loc_Hz)
            F_idx_low = max(0, F_idx - 2)
            F_idx_high = min(len(aid.fitter.F_Hz) - 1, F_idx + 2)
            #get the local mean width of bins to restrict bandwidths
            F_mean_Hz = (aid.fitter.F_Hz[F_idx_high] - aid.fitter.F_Hz[F_idx_low]) / (F_idx_high - F_idx_low)
            if b.Q_rank_BW < F_mean_Hz / 2:
                exclude_p.add(b.idx_p)
                exclude_z.add(b.idx_z)
                duals_used.append(b)
                aid.log(F_mean_Hz, b.Q_rank_BW)
                aid.log("Removing by bandwidth")
                used = True

        if not used:
            duals_rejected.append(b)

    #switch method to a residual ranking
    duals.sort(key = lambda b: b.F_loc_Hz)
    if False and not duals_used:
        duals_rejected = []
        for b in duals:
            used = False
            add_resrank(aid.fitter, b)
            if b.idx_max - b.idx_min < 2:
                #kill single point pairs (they shouldn't exist anyway)
                exclude_p.add(b.idx_p)
                exclude_z.add(b.idx_z)
                duals_used.append(b)
                used = True

            if used:
                pass
            if not used:
                duals_rejected.append(b)

    for idx_p, coding in enumerate(aid.fitter.den_codings):
        if idx_p not in exclude_p:
            den_codings.append(coding)
    for idx_z, coding in enumerate(aid.fitter.num_codings):
        if idx_z not in exclude_z:
            num_codings.append(coding)

    fitter_out = aid.fitter.__class__(
        parent            = aid.fitter,
        num_codings       = num_codings,
        den_codings       = den_codings,
    )

    return declarative.Bunch(locals())


def order_reduce_weakroots(
    aid,
    delay_s,
):
    """
    """

    aid = ensure_aid(aid)

    F_largest_Hz = np.max(aid.fitter.F_Hz)
    Xsf = 1j * F_largest_Hz
    def ang_mag_mapping(coding_list, rtype):
        new_list = []
        for coding in coding_list:
            try:
                roots = coding.roots_Sf()
            except fitters_ZPK.BranchCutAmbiguity:
                continue
            if np.any(abs(np.real(roots)) < 1 * F_largest_Hz):
                new_list.append(None)
                continue
            #aid.log(roots, rtype)
            #divide by the distance to the roots to normalize DC to 1
            TF = TFmath.polyvalfromroots(Xsf, roots) / np.prod(np.abs(roots))
            mag = abs(TF)
            angle = TF / mag * np.prod(-np.sign(roots))
            #print(roots, mag, angle)
            new_list.append((angle, np.log(mag)))
        #the LAST entry is the null-map, so it must be checked-for
        new_list.append((1, 0))
        return new_list
    z_map = ang_mag_mapping(aid.fitter.num_codings, 'zeros')
    p_map = ang_mag_mapping(aid.fitter.den_codings, 'poles')
    #print(z_map, p_map)

    def nearest_idx(lst_1, lst_2, num_first = True):
        nearest_lst = []
        for idx_1, pair_1 in enumerate(lst_1):
            if pair_1 is None:
                continue
            ang_1, mag_1 = pair_1
            dist_nearest = float('inf')
            ang_diff_nearest = float('inf')
            mag_diff_nearest = float('inf')
            idx_nearest = None
            for idx_2, pair_2 in enumerate(lst_2):
                if pair_2 is None:
                    continue
                if (idx_2 == len(lst_2) - 1) and (idx_1 == len(lst_1) - 1):
                    continue
                ang_2, mag_2 = pair_2
                ang_diff = np.angle(ang_1 * ang_2.conjugate())
                mag_diff = mag_1 - mag_2
                dist = (ang_diff**2 + 1 * mag_diff**2)**.5
                if dist < dist_nearest:
                    idx_nearest = idx_2
                    dist_nearest = dist
                    ang_diff_nearest = ang_diff
                    mag_diff_nearest = mag_diff
            if idx_nearest is None:
                continue
            if idx_1 < len(lst_1) - 1:
                if num_first:
                    coding_num = aid.fitter.num_codings[idx_1]
                    idx_num = idx_1
                else:
                    coding_den = aid.fitter.den_codings[idx_1]
                    idx_den = idx_1
            else:
                if num_first:
                    coding_num = None
                    idx_num = None
                else:
                    coding_den = None
                    idx_den = None
            if idx_nearest < len(lst_2) - 1:
                if not num_first:
                    coding_num = aid.fitter.num_codings[idx_nearest]
                    idx_num = idx_nearest
                else:
                    coding_den = aid.fitter.den_codings[idx_nearest]
                    idx_den = idx_nearest
            else:
                if not num_first:
                    coding_num = None
                    idx_num = None
                else:
                    coding_den = None
                    idx_den = None
            data = declarative.Bunch(
                idx_nearest = idx_nearest,
                coding_num = coding_num,
                coding_den = coding_den,
                idx_num = idx_num,
                idx_den = idx_den,
                dist = dist_nearest,
                ang_diff = ang_diff_nearest,
                mag_diff = mag_diff_nearest,
            )
            nearest_lst.append(data)
        return nearest_lst
    z_nearest = nearest_idx(z_map, p_map, num_first = True)
    p_nearest = nearest_idx(p_map, z_map, num_first = False)

    key = lambda b : b.dist if b is not None else float('inf')
    if z_nearest:
        z_min = min(z_nearest, key = key)
    else:
        z_min = None
    if p_nearest:
        p_min = min(p_nearest, key = key)
    else:
        p_min = None
    if p_min is None:
        if z_min is None:
            return False
        else:
            zp_min = z_min
    else:
        if z_min is None:
            zp_min = p_min
        else:
            zp_min = min(z_min, p_min, key = key)

    num_codings = list(aid.fitter.num_codings)
    den_codings = list(aid.fitter.den_codings)
    if zp_min.idx_num is not None:
        del num_codings[zp_min.idx_num]
    if zp_min.idx_den is not None:
        del den_codings[zp_min.idx_den]

    if zp_min.coding_num is not None:
        pr_num = zp_min.coding_num.roots_Sf()
    else:
        pr_num = None

    if zp_min.coding_den is not None:
        pr_den = zp_min.coding_den.roots_Sf()
    else:
        pr_den = None

    aid.log("WEAK REMOVE: ", pr_num, pr_den)
    fitter_new = aid.fitter.copy(
        num_codings = num_codings,
        den_codings = den_codings,
    )
    ratio = fitter_new.residuals_average / aid.fitter.residuals_average
    aid.log("RATIO for WEAK1: ", ratio)
    #aid.log("BEFORE fitters_ZPK: ", aid.fitter.residuals_average)
    #aid.log("BEFORE NEW: ", fitter_new.residuals_average)

    #maybe feed forward the delay
    #fitter_new.optimize_delay(only_delay = False)
    if delay_s is None:
        fitter_new.optimize_delay(
            residuals_type = aid.hint('residuals_type_alt', default = 'log'),
            only_delay = False,
            rescale = True,
            aid = aid,
        )
        fitter_new.optimize_delay(
            rescale = True,
            only_delay = False,
            aid = aid,
        )
    else:
        fitter_new.optimize(
            residuals_type = aid.hint('residuals_type_alt', default = 'log'),
            rescale = True,
            aid = aid,
        )
        fitter_new.optimize(
            rescale = True,
            aid = aid,
        )

    ratio = fitter_new.residuals_average / aid.fitter.residuals_average
    aid.log("RATIO for WEAK2: ", ratio)
    algorithms.sign_check_flip(fitter_new)
    return aid.fitter_check(
        fitter_new,
        hint_name = 'weakroots',
        variant = 'OrdDn',
    )


def order_reduce_aggressive(aid):
    aid = ensure_aid(aid)

    N = 1
    #DOESN'T SEEM TO HELP
    #duals = match_pairs(aid.fitter)
    #ignore_orig = aid.fitter.codings_ignore
    #for dualb in duals:
    #    aid.fitter.codings_ignore = aid.fitter.root_codings_set - set([
    #        aid.fitter.num_codings[dualb.idx_z],
    #        aid.fitter.den_codings[dualb.idx_p]
    #    ]).union(ignore_orig)
    #    aid.fitter.optimize()
    #aid.fitter.codings_ignore = ignore_orig
    while True:
        duals = match_pairs(aid = aid)
        #TODO annotate better
        #at this cutoff, the pairs are essentially guaranteed not to affect the residual in a meaningful way
        ORb = order_reduce(
            aid = aid,
            duals = duals,
            Q_rank_cutoff = .7,
        )
        def closure(ORb):
            def plot_duals(fitter):
                axB = plots.plot_fit(fitter)
                plots.plot_ZP_grab(
                    fitter,
                    ORb.duals_used,
                    axB = axB,
                    color = 'black'
                )
                plots.plot_ZP_grab(
                    fitter,
                    ORb.duals_rejected,
                    axB = axB,
                    color = 'blue'
                )
                return axB
            return plot_duals
        aid.annotate(
            name = 'duals_ID_{0}'.format(N),
            method = 'MultiReprFilter.optimize',
            fitter   = aid.fitter,
            plotter  = closure(ORb),
            verbosity = 4,
            about = (
                """
                ID Pairs for order reduction
                """
            ),
        )
        if ORb.duals_used:
            aid.fitter_update(ORb.fitter_out)

        further_improvement = False
        if not ORb.duals_used:
            duals = match_pairs(duals = False, aid = aid)
            duals = [add_resrank(b, aid = aid) for b in duals]
            duals.sort(key = lambda b : b.res_h_tot)
            for b_low in duals[:4]:
                #aid.log(b_low.res_h_tot, b_low.res_mod_tot, b_low.res_tot)
                num_codings = list(aid.fitter.num_codings)
                den_codings = list(aid.fitter.den_codings)
                #aid.log("LEN:", len(aid.fitter.num_codings), len(aid.fitter.den_codings))
                #aid.log("IDX:", b_low.idx_z, b_low.idx_p)
                del num_codings[b_low.idx_z]
                del den_codings[b_low.idx_p]

                fitter_new = aid.fitter.__class__(
                    parent = aid.fitter,
                    num_codings = num_codings,
                    den_codings = den_codings,
                )
                fitter_new.optimize(aid = aid)
                res_ratio = fitter_new.residuals_average / aid.fitter.residuals_average
                aid.log("RATIO: ", res_ratio)
                if res_ratio < 1.0:
                    aid.log("fit improved from pair at", b_low.F_loc_Hz)
                    further_improvement = True
                    aid.fitter_update(fitter_new)
                    break
                else:
                    aid.log("fit NOT improved from pair at", b_low.F_loc_Hz)
        if not ORb.duals_used and not further_improvement:
            break

        aid.fitter.optimize(aid = aid)
        #annotate.annotate(
        #    doc_db,
        #    name = 'seq_order_reduce_opt_{0}'.format(N),
        #    method = 'MultiReprFilter.optimize',
        #    fitter   = aid.fitter,
        #    plotter  = plots.plot_fit,
        #    verbosity = 6,
        #    about = (
        #        """
        #        Optimize again after order reduction
        #        """
        #    ),
        #)
        N += 1
        aid.log("N: ", N)
    return


def order_reduce_negroots(aid):
    """
    ONLY works on F_nyquist_Hz > 2 * F_max_Hz
    """
    aid = ensure_aid(aid)

    #Only works in Z domain
    assert(aid.fitter.F_nyquist_Hz is not None)

    #F_largest_Hz = np.max(aid.fitter.F_Hz)
    def BW_check(coding_list, rtype):
        new_list = []
        for coding in coding_list:
            root_list = np.asarray(coding.roots())
            #exclude positive ones
            if np.any(root_list.real > 0):
                new_list.append(coding)
                continue
            aid.log("REMOVING ", rtype, root_list)
        return new_list

    num_codings = BW_check(aid.fitter.num_codings, 'Z: ')
    den_codings = BW_check(aid.fitter.den_codings, 'P: ')

    fitter_new = aid.fitter.__class__(
        parent      = aid.fitter,
        num_codings = num_codings,
        den_codings = den_codings,
    )
    aid.fitter_update(fitter_new)
    return
