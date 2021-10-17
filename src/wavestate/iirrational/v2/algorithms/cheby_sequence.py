# -*- coding: utf-8 -*-
"""
"""


from ...fitters_rational import ChebychevFilter
from ...utilities import ensure_aid


def reldeg_split(relative_degree):
    if relative_degree is None or relative_degree == 0:
        reln = 0
        reld = 0
    elif relative_degree > 0:
        reln = relative_degree
        reld = 0
    else:
        reln = 0
        reld = -relative_degree
    return reln, reld


def rational_cheby_fit(
    ZPKrep,
    order           = None,
    order_max       = None,
    order_min       = 20,
    relative_degree = None,
    aid             = None,
):
    aid              = ensure_aid(aid)

    if ZPKrep.F_nyquist_Hz is not None:
        raise RuntimeError("v2 only supports S domain for now")

    if order is not None:
        order     = int(order)
    order_max = int(order_max)
    order_min = int(order_min)

    N_first = order_min
    N_final = int(min(len(ZPKrep.F_Hz) // 10, order_max))

    if order is not None:
        return cheby_single(
            ZPKrep,
            relative_degree  = relative_degree,
            order            = order,
            aid              = aid,
        )

    #otherwise, scan through
    N_current = N_first

    fitter_last = cheby_single(
        ZPKrep,
        order            = N_current ,
        relative_degree  = relative_degree,
        aid              = aid,
    )
    restot_last = fitter_last.residuals_average

    def fitter_ord(fitter):
        return max(len(fitter.zeros), len(fitter.poles))

    with aid.log_heading("Rational Fit Order Scanning"):
        while True:
            if N_current == N_final:
                fitter_use = fitter_last
                break

            N_current = N_current * 2
            if N_current > N_final:
                N_current = N_final

            fitter = cheby_single(
                ZPKrep,
                order           = N_current,
                relative_degree = relative_degree,
                aid             = aid,
            )
            restot = fitter.residuals_average
            fitter_red = fitter.copy()
            fitter_red.matched_pairs_clear(Q_rank_cutoff = .2)
            restot_red = fitter_red.residuals_average

            if restot_last < restot:
                fitter_use = fitter_last
                aid.log_debug(8, "Using last (direct)!", fitter_last.order)
                break

            if restot_last < restot_red:
                fitter_use = fitter_last
                aid.log_debug(8, "Using last (reduced)!", fitter_last.order)
                break

            if restot_last < 1.10 * restot:
                ord_red = fitter_red.order
                ord_last = fitter_last.order
                if ord_red < ord_last:
                    fitter_use = fitter_red
                else:
                    fitter_use = fitter_last
                aid.log_debug(8, "Using current")
                break
            #else continue the loop
            fitter_last = fitter
            restot_last = restot

    return fitter_use


def cheby_single(
    ZPKrep,
    order           = 20,
    relative_degree = 0,
    aid             = None,
):
    aid = ensure_aid(aid)

    if ZPKrep.F_nyquist_Hz is not None:
        raise RuntimeError("F_nyquist_Hz must be None (ZPKrep must be in S domain)")

    if relative_degree is None:
        relative_degree = 0

    reln, reld = reldeg_split(relative_degree)
    fitter = ChebychevFilter(
        ZPKrep = ZPKrep,
        nzeros = (order//2) + reln,
        npoles = (order//2) + reld,
    )

    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.nzeros = order + reln
    fitter.npoles = order + reld
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
    #fitter.stabilize = False
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    return fitter


