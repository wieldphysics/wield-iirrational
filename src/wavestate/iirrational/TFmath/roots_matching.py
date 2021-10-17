# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import declarative


def nearest_idx(
        lst_1,
        lst_2 = None,
        metric_pair_dist = None,
        return_distances = False,
):
    """
    If lst_2 is given, this returns all of the nearest items in lst_2 to lst_1.
    If not given, this returns all of the nearest elements of lst_1 to itself,
    ignoring self elements.

    if metric_pair_dist is None, use the standard distance on complex plane.
    This is the fastest.
    """
    dists = []
    if lst_2 is not None:
        #TODO, this could be much more efficient with sorting..
        if metric_pair_dist is None:
            def metric_pair_dist(r1, r2):
                return abs(r1 - r2)
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
                dist = metric_pair_dist(r1, r2)
                if dist < dist_nearest:
                    idx_nearest = idx_2
                    dist_nearest = dist
            nearest_lst.append(idx_nearest)
            dists.append(dist_nearest)
    else:
        #TODO, this could be much more efficient with sorting..
        if metric_pair_dist is None:
            def metric_pair_dist(r1, r2):
                return abs(r1 - r2)
        nearest_lst = []
        for idx_1, r1 in enumerate(lst_1):
            if r1 is None:
                nearest_lst.append(None)
                continue
            dist_nearest = float('inf')
            idx_nearest = None
            for idx_2, r2 in enumerate(lst_1):
                if idx_2 == idx_1:
                    continue
                if r2 is None:
                    continue
                dist = metric_pair_dist(r1, r2)
                if dist < dist_nearest:
                    idx_nearest = idx_2
                    dist_nearest = dist
            nearest_lst.append(idx_nearest)
            dists.append(dist_nearest)
    if return_distances:
        return nearest_lst, dists
    else:
        return nearest_lst


def nearest_unique_idx(l1, l2):
    idx_p_used = []
    idx_list = []
    l1 = list(l1)

    for idx_1, idx_2 in enumerate(l1):
        if idx_2 is None:
            continue
        if idx_1 == l2[idx_2]:
            idx_p_used.append(idx_2)
            idx_list.append(
                (idx_1, idx_2)
            )
        else:
            l1[idx_1] = None
    l2_copy = [None] * len(l2)
    for idx_2 in idx_p_used:
        l2_copy[idx_2] = l2[idx_2]
    return wavestate.bunch.Bunch(
        idx_list = idx_list,
        l1       = l1,
        l2       = l2_copy,
    )


def nearest_unique_pairs(
    l1,
    l2,
    metric_pair_dist = None,
):
    r12_list = []
    idx_list  = []
    l1 = list(l1)
    l2 = list(l2)
    l1_nearest = nearest_idx(l1, l2, metric_pair_dist = metric_pair_dist)
    l2_nearest = nearest_idx(l2, l1, metric_pair_dist = metric_pair_dist)

    l1_remain = []
    l2_remain = []
    idx_2_used = []

    for idx_1, idx_2 in enumerate(l1_nearest):
        if idx_2 is None:
            l1_remain.append(l1[idx_1])
            continue
        #coding_z = aid.fitter.num_codings[idx_1]
        #coding_p = aid.fitter.den_codings[idx_2]
        #TODO annotate about stability
        p = l2[idx_2]
        z = l1[idx_1]
        if idx_1 == l2_nearest[idx_2]:
            idx_2_used.append(idx_2)
            r12_list.append(
                (z, p)
            )
            idx_list.append(
                (idx_1, idx_2)
            )
        else:
            l1_remain.append(l1[idx_1])
            l1_nearest[idx_1] = None
    idx_2_used = set(idx_2_used)
    for idx_2, p in enumerate(l2):
        if idx_2 not in idx_2_used:
            l2_remain.append(p)
            l2_nearest[idx_2] = None
    assert(len(r12_list) + len(l1_remain) == len(l1))
    assert(len(r12_list) + len(l2_remain) == len(l2))
    return wavestate.bunch.Bunch(
        r12_list   = r12_list,
        l1_remain  = l1_remain,
        l2_remain  = l2_remain,
        idx_list   = idx_list,
        l1         = l1_nearest,
        l2         = l2_nearest,
    )

def nearest_pairs(
    l1,
    l2,
    metric_pair_dist = None,
):
    #TODO, allow other rankings than distance

    rpB = nearest_unique_pairs(l1, l2, metric_pair_dist)
    #not going to maintain these lists
    del rpB.idx_list
    del rpB.l1
    del rpB.l2

    while True:
        pair_lists = []
        l1_nearest, l1_dist = nearest_idx(
            rpB.l1_remain,
            rpB.l2_remain,
            metric_pair_dist = metric_pair_dist,
            return_distances = True,
        )
        for idx_1, idx_2 in enumerate(l1_nearest):
            if idx_2 is None:
                continue
            dist = l1_dist[idx_1]
            pair_lists.append((dist, idx_1, idx_2))
        l2_nearest, l2_dist = nearest_idx(
            rpB.l2_remain,
            rpB.l1_remain,
            metric_pair_dist = metric_pair_dist,
            return_distances = True,
        )
        for idx_2, idx_1 in enumerate(l2_nearest):
            if idx_1 is None:
                continue
            dist = l2_dist[idx_2]
            pair_lists.append((dist, idx_1, idx_2))
        if not pair_lists:
            break
        pair_lists.sort()
        dist, idx_1, idx_2 = pair_lists[0]
        rpB.r12_list.append((rpB.l1_remain[idx_1], rpB.l2_remain[idx_2]))
        del rpB.l1_remain[idx_1]
        del rpB.l2_remain[idx_2]
    return rpB

def SOS_pair_rolloff(
    Zr,
    Zc,
    Pr,
    Pc,
    metric_rolloff = abs,
):
    metric = metric_rolloff
    zzpp_list = []

    Zr = list(Zr)
    Zc = list(Zc)
    Pr = list(Pr)
    Pc = list(Pc)

    #sort on their rolloff bandwidth
    Zc.sort(key = metric)
    Zr.sort(key = metric)
    Pc.sort(key = metric)
    Pr.sort(key = metric)

    def real_low(Rr, Rc):
        if Rc:
            if Rr:
                if metric(Rr[0]) < metric(Rc[0]):
                    r_low = True
                else:
                    r_low = False
            else:
                r_low = False
        else:
            r_low = True
        return r_low

    def real_hold_insert(rpair, rhold):
        if rhold is None:
            return rpair
        z1, p1 = rhold
        z2, p2 = rpair
        zzpp_list.append((z1, z2, p1, p2))
        return None

    rpair_hold = None
    while True:
        Pr_low = real_low(Pr, Pc)
        Zr_low = real_low(Zr, Zc)

        if Pr_low and Zr_low:
            if Pr and Zr:
                p = Pr.pop(0)
                z = Zr.pop(0)
                rpair_hold = real_hold_insert((z, p), rpair_hold)
                continue
            elif Pr:
                #can only occur if Zc is also empty
                p = Pr.pop(0)
                rpair_hold = real_hold_insert((None, p), rpair_hold)
                continue
            elif Zr:
                #can only occur if Pc is also empty
                z = Zr.pop(0)
                rpair_hold = real_hold_insert((z, None), rpair_hold)
                continue
            else:
                #can only occur when all lists are empty
                break

        elif not Pr_low and not Zr_low:
            if Pc and Zc:
                p = Pc.pop(0)
                z = Zc.pop(0)
                zzpp_list.append((z, z.conjugate(), p, p.conjugate()))
                continue
            elif Pc:
                #can only occur if Zr is also empty
                p = Pc.pop(0)
                zzpp_list.append((None, None, p, p.conjugate()))
                continue
            elif Zc:
                #can only occur if Pr is also empty
                z = Zc.pop(0)
                zzpp_list.append((z, z.conjugate(), None, None))
                continue
            else:
                #can only occur when all lists are empty
                break

        #now it must be a mix for Pr_low and Zr_low, so work around the complex one
        elif Pr_low:
            #must be that Zr_low is False, Pr_low is true
            z = Zc.pop(0)
            if Pc:
                pc_t = Pc[0]
                if len(Pr) > 1:
                    pr_t = Pr[1]
                    if metric(pr_t) < metric(pc_t):
                        #we have a DOUBLE low real, so use that
                        p1 = Pr.pop(0)
                        p2 = Pr.pop(0)
                        zzpp_list.append((z, z.conjugate(), p1, p2))
                        continue
                    else:
                        #complex is between the others, so use it instead
                        p = Pc.pop(0)
                        zzpp_list.append((z, z.conjugate(), p, p.conjugate()))
                        continue
                else:
                    #if Pr is almost empty, then favor the complex
                    p = Pc.pop(0)
                    zzpp_list.append((z, z.conjugate(), p, p.conjugate()))
                    continue
            else:
                #no more complex poles, so must exhaust real ones
                if len(Pr) > 1:
                    p1 = Pr.pop(0)
                    p2 = Pr.pop(0)
                    zzpp_list.append((z, z.conjugate(), p1, p2))
                    continue
                elif len(Pr) == 1:
                    p1 = Pr.pop(0)
                    zzpp_list.append((z, z.conjugate(), p1, None))
                    continue
                else:
                    zzpp_list.append((z, z.conjugate(), None, None))
                    continue
        else:
            #must be that Pr_low is False
            #must be that Zr_low is False, Pr_low is true
            p = Pc.pop(0)
            if Zc:
                zc_t = Zc[0]
                if len(Zr) > 1:
                    zr_t = Zr[1]
                    if metric(zr_t) < metric(zc_t):
                        #we have a DOUBLE low real, so use that
                        z1 = Zr.pop(0)
                        z2 = Zr.pop(0)
                        zzpp_list.append((z1, z2, p, p.conjugate()))
                        continue
                    else:
                        #complex is between the others, so use it instead
                        z = Zc.pop(0)
                        zzpp_list.append((z, z.conjugate(), p, p.conjugate()))
                        continue
                else:
                    #if Zr is almost empty, then favor the complex
                    z = Zc.pop(0)
                    zzpp_list.append((z, z.conjugate(), p, p.conjugate()))
                    continue
            else:
                #no more complex zeros, so must exhaust real ones
                if len(Zr) > 1:
                    z1 = Zr.pop(0)
                    z2 = Zr.pop(0)
                    zzpp_list.append((z1, z2, p, p.conjugate()))
                    continue
                elif len(Zr) == 1:
                    z1 = Zr.pop(0)
                    zzpp_list.append((z1, None, p, p.conjugate()))
                    continue
                else:
                    zzpp_list.append((None, None, p, p.conjugate()))
                    continue

    if rpair_hold is not None:
        z, p = rpair_hold
        zzpp_list.append((z, None, p, None))

    return zzpp_list


def match_SOS_pairs(
    Zr,
    Zc,
    Pr,
    Pc,
    F_nyquist_Hz = None,
    metric_rolloff = None,
    metric_pair_dist = None,
):
    """
    Match and create pairs suitable for SOS representation. The output
    is a list of 4-tuples with z1, z2, p1, p2. If roots are complex,
    they are guaranteed to be partnered with their conjugate.

    Some z, or p may be None, indicating that the system ran out.
    """

    Nz_tot = len(Zc) * 2 + len(Zr)
    Np_tot = len(Pc) * 2 + len(Pr)

    if metric_pair_dist is None:
        if F_nyquist_Hz is None:
            def metric_pair_dist(r1, r2):
                return abs(r1 - r2)
        else:
            def metric_pair_dist(r1, r2):
                return abs(r1 - r2)

    pairB = nearest_unique_pairs(Zc, Pc, metric_pair_dist = metric_pair_dist)
    Zc = pairB.l1_remain
    Pc = pairB.l2_remain

    zzpp_list = []
    for z, p in pairB.r12_list:
        zzpp_list.append(
            (z, z.conjugate(), p, p.conjugate())
        )

    if metric_rolloff is None:
        if Np_tot < Nz_tot:
            #when the filter is AC coupled, apply strongest
            #AC coupling early and weaker later
            if F_nyquist_Hz is None:
                def metric_rolloff(r):
                    return abs(r)
            else:
                def metric_rolloff(r):
                    F = np.angle(r)
                    d = 1 - abs(r)
                    return (d**2 + F**2)**.5
        else:
            #when the filter has rolloff, apply the strongest rolloff later
            if F_nyquist_Hz is None:
                def metric_rolloff(r):
                    return 1/abs(r)
            else:
                def metric_rolloff(r):
                    F = np.angle(r)
                    d = 1 - abs(r)
                    return 1 / (d**2 + F**2)**.5

    zzpp_list.extend(
        SOS_pair_rolloff(Zr, Zc, Pr, Pc, metric_rolloff = metric_rolloff)
    )

    Nz = 0
    Np = 0
    for (z1, z2, p1, p2) in zzpp_list:
        if z1 is not None:
            Nz += 1
        if z2 is not None:
            Nz += 1
        if p1 is not None:
            Np += 1
        if p2 is not None:
            Np += 1
    assert(Nz == Nz_tot)
    assert(Np == Np_tot)
    return zzpp_list

