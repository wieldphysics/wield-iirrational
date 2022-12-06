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


def fit_AAA(aid, order_hint=None):
    if order_hint is None:
        return fit_AAA_base(
            aid,
            order=aid.hint(
                "rational_AAA_fit_order",
                "rational_fit_order",
                "order_initial",
            ),
            order_max=aid.hint(
                "rational_AAA_fit_order_max",
                "rational_fit_order_max",
                "order_max",
            ),
            order_min=aid.hint(
                "rational_AAA_fit_order_min",
                "rational_fit_order_min",
                "order_min",
            ),
        )
    else:
        return fit_AAA_base(
            aid,
            order=aid.hint(
                order_hint,
                "rational_AAA_fit_order",
                "rational_fit_order",
                "order_initial",
            ),
            order_max=aid.hint(
                "rational_AAA_fit_order_max",
                "rational_fit_order_max",
                "order_max",
            ),
            order_min=aid.hint(
                "rational_AAA_fit_order_min",
                "rational_fit_order_min",
                "order_min",
            ),
        )


def fit_AAA_base(
    aid,
    order,
    order_max,
    order_min,
    account_size=True,
):
    if order == 0:
        return

    aid.log_progress(4, "AAA rational fit")
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

    from wavestate.control.AAA import tfAAA
    if order is None:
        order = 20

    if order_max is None:
        order_max = order
    else:
        print(order, order_max)
        order_max = min(order, order_max)

    if account_size:
        factors_order = aid.fitter_orders().factors_maxzp
        # print("FACTORS ORDER", factors_order)
        order_max = order_max - factors_order
        order_max = max(order_max, 6)

    aaa = tfAAA(
        aid.fitter.F_Hz,
        aid.fitter.data_no_overlay,
        exact=False,
        res_tol=None,
        s_tol=None,
        w=aid.fitter.W,
        w_res=None,
        degree_max=order_max,
        nconv=None,
        nrel=10,
        rtype="log",
        lf_eager=True,
        supports=(),
        minreal_cutoff=None,
    )
    zeros = aaa.zeros# / (2 * np.pi)
    poles = aaa.poles# / (2 * np.pi)
    select = poles.real > 0
    poles[select] = -poles[select].conjugate()
    gain = aaa.gain# * (2 * np.pi)**(len(zeros) - len(poles))
    gain = (-1)**np.sum(select) * gain
    fitter = aid.fitter.regenerate(
        ZPKrep=aid.fitter.ZPKrep,
        zeros=zeros,
        poles=poles,
        gain=gain
    )
    print(np.all(fitter.poles.fullplane.real < 0))

    # TODO, add debug_AAA hint
    # from .. import plots
    # axB = plots.plot_fitter_flag(fitter=aid.fitter, xscale='log')
    # axB.save("AAA_pre_{}.pdf".format(aid.N_update))
    # axB = plots.plot_fitter_flag(fitter=fitter, xscale='log')
    # axB.save("AAA_{}.pdf".format(aid.N_update))

    # print("AAAlist: ", poles, zeros)

    aid.fitter_update(
        fitter,
        representative=False,
    )
    return

