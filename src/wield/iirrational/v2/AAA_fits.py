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

    from wield.control.AAA import tfAAA
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

    supports = np.linspace(aid.fitter.F_Hz[0], aid.fitter.F_Hz[-1], 5)
    rsum = np.cumsum(aid.fitter.residuals_sq)
    rsum = rsum / rsum[-1]
    sratios = np.linspace(0, 1, 8)[:-1]
    idx = np.searchsorted(rsum, sratios)
    supports = aid.fitter.F_Hz[idx]

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
        supports=supports,
        minreal_cutoff=None,
        all_real=False
    )


    from wield.control import SISO

    #order = aaa.order
    #order_orig = order
    #while order > 1:
    #    aaa.choose(order)
    #    zeros = aaa.zeros
    #    poles = aaa.poles
    #    gain = aaa.gain
    #    select = poles.real > 0
    #    num_unstable = np.sum(select)
    #    
    #    aid.log_progress(4, "AAA order {}, num unstable: {} of {}, residuals {:.2f}".format(order, num_unstable, len(poles), aaa.fit_dict['res_rms']**2))
    #    if num_unstable == 0:
    #        break
    #    order -= 1
    #if order == 0:
    #    aaa.choose(order_orig)
    #    aid.log_progress(4, "AAA always unstable, using order {}".format(aaa.order))

    zeros = aaa.zeros
    poles = aaa.poles
    gain = aaa.gain

    from .. import fitters_ZPK
    fitter_bad = aid.fitter.regenerate(
        ZPKrep=aid.fitter.ZPKrep,
        coding_map=fitters_ZPK.codings_s.coding_maps.SOS,
        poles=np.asarray(poles),
        zeros=np.asarray(zeros),
        gain=gain,
        check_sign=True,
    )


    #from wield.utilities.mpl import mplfigB
    #axB = mplfigB(Nrows=2)
    #
    #h = fitter_bad.data_no_overlay
    #axB.ax0.loglog(fitter_bad.F_Hz, abs(h), marker='.', ls='')
    #axB.ax1.semilogx(fitter_bad.F_Hz, np.angle(h, deg=True), marker='.', ls='')

    #h = aaa(fitter_bad.F_Hz)
    #axB.ax0.loglog(fitter_bad.F_Hz, abs(h))
    #axB.ax1.semilogx(fitter_bad.F_Hz, np.angle(h, deg=True))


    #h = fitter_bad.xfer_fit_no_overlay
    #axB.ax0.loglog(fitter_bad.F_Hz, abs(h))
    #axB.ax1.semilogx(fitter_bad.F_Hz, np.angle(h, deg=True))

    # TODO, add debug_AAA hint
    # from .. import plots
    # axB = plots.plot_fitter_flag_residuals(fitter=aid.fitter, xscale='log')
    # axB.save("AAA_pre_{}.pdf".format(aid.N_update))
    # axB = plots.plot_fitter_flag_residuals(fitter=fitter_bad, xscale='log')
    # axB.save("AAA_{}.pdf".format(aid.N_update))

    with fitter_bad.with_codings_only([fitter_bad.gain_coding]):
        fitter_bad.optimize()
    fitter_bad.optimize()

    poles = fitter_bad.poles.fullplane
    select = poles.real > 0
    poles[select] = -poles[select].conjugate()
    gain = (-1)**np.sum(select) * gain

    fitter = aid.fitter.regenerate(
        ZPKrep=aid.fitter.ZPKrep,
        zeros=fitter_bad.zeros,
        poles=poles,
        gain=fitter_bad.gain,
    )
    fitter.optimize()
    assert(np.all(fitter.poles.fullplane.real < 0))

    #h = fitter_bad.xfer_fit_no_overlay
    #axB.ax0.loglog(fitter_bad.F_Hz, abs(h))
    #axB.ax1.semilogx(fitter_bad.F_Hz, np.angle(h, deg=True))

    #h = fitter.xfer_fit_no_overlay
    #axB.ax0.loglog(fitter.F_Hz, abs(h))
    #axB.ax1.semilogx(fitter.F_Hz, np.angle(h, deg=True))

    #for z in aaa.supports:
    #    axB.ax0.axvline(z, lw=0.5, ls='--', color='black')

    #axB.save("AAA_tests_{}.pdf".format(aid.N_update))

    #print("AAAlist: ", poles, zeros)

    aid.fitter_update(
        fitter,
        representative=False,
    )
    return

