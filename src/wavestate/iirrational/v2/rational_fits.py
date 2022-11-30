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

from .algorithms import (
    cheby_sequence,
    disc_sequence,
)



def fit_cheby(aid, order_hint=None):
    if order_hint is None:
        return fit_cheby_base(
            aid,
            order=aid.hint(
                "rational_cheby_fit_order",
                "rational_fit_order",
                "order_initial",
            ),
            order_max=aid.hint(
                "rational_cheby_fit_order_max",
                "rational_fit_order_max",
                "order_max",
            ),
            order_min=aid.hint(
                "rational_cheby_fit_order_min",
                "rational_fit_order_min",
                "order_min",
            ),
        )
    else:
        return fit_cheby_base(
            aid,
            order=aid.hint(
                order_hint,
                "rational_cheby_fit_order",
                "rational_fit_order",
                "order_initial",
            ),
            order_max=aid.hint(
                "rational_cheby_fit_order_max",
                "rational_fit_order_max",
                "order_max",
            ),
            order_min=aid.hint(
                "rational_cheby_fit_order_min",
                "rational_fit_order_min",
                "order_min",
            ),
        )


def fit_cheby_base(
    aid,
    order,
    order_max,
    order_min,
):
    if order == 0:
        return

    aid.log_progress(4, "chebychev rational fit")
    # rat_fitter = v2.disc_sequence.rational_disc_fit(
    reldeg = aid.hint("relative_degree")
    if reldeg is None:
        reldeg_max = aid.hint("relative_degree_max")
        reldeg_min = aid.hint("relative_degree_min")
        if reldeg_min is None and reldeg_max is None:
            reldeg = None
        elif reldeg_min is None:
            reldeg = 0
        elif reldeg_max is None:
            reldeg = 0
        else:
            reldeg = int((reldeg_min + reldeg_max) // 2)

    factor_orders = aid.fitter_orders("factors")
    if reldeg is not None:
        diff_reldeg = reldeg - factor_orders.reldeg
    else:
        diff_reldeg = 0

    factors_order = aid.fitter_orders().factors_maxzp
    print("RELDEG: ", factor_orders.reldeg, aid.fitter_orders().factors_reldeg, diff_reldeg)
    # print("FACTORS ORDER", factors_order)
    if order is not None:
        order -= factors_order
        order = max(order, 4)

    rat_fitter = cheby_sequence.rational_cheby_fit(
        aid=aid,
        order_max=order_max - factors_order,
        order_min=max(order_min - factors_order, 2),
        order=order,
        relative_degree=diff_reldeg,
        ZPKrep=aid.fitter.ZPKrep,
    )

    fitter = aid.fitter.regenerate(
        ZPKrep=rat_fitter.ZPKrep,
    )

    aid.fitter_update(
        fitter,
        representative=True,
    )


def fit_disc_factored(aid):
    with aid.factorization():
        fit_disc(aid)


def fit_disc(aid):
    aid.log_progress(4, "rational disc fitting")
    order = aid.hint(
        "rational_cheby_fit_order",
        "rational_fit_order",
        "order_initial",
    )
    order_max = aid.hint(
        "rational_cheby_fit_order_max",
        "rational_fit_order_max",
        "order_max",
    )
    order_min = aid.hint(
        "rational_cheby_fit_order_min",
        "rational_fit_order_min",
        "order_min",
    )

    # TODO
    from wavestate.iirrational.v1 import disc_sequence

    rat_fitter = disc_sequence.rational_disc_fit(
        F_Hz=aid.fitter.F_Hz,
        data=aid.fitter.data,
        nyquist_final_Hz=None,
        SNR=aid.fitter.W,
        order=order,
        order_max=order_max,
    )
    # rat_fitter = disc_sequence.rational_disc_fit(
    #    aid             = aid,
    #    order_max       = order_max,
    #    order_min       = order_min,
    #    order           = order,
    #    ZPKrep          = aid.fitter.ZPKrep,
    # )

    fitter = aid.fitter.regenerate(
        ZPKrep=rat_fitter.ZPKrep,
    )

    # TODO
    # validate = aid.hint('rational_fitter_validate')
    # if validate is not None:
    #    validate(rat_fitter, fitter)

    aid.fitter_update(
        fitter,
        representative=True,
    )

    aid.log(
        1,
        "Initial Order: (Z={0}, P={1}, Z-P={2})".format(
            len(fitter.zeros),
            len(fitter.poles),
            len(fitter.zeros) - len(fitter.poles),
        ),
    )
