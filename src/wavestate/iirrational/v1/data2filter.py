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
from wavestate import declarative

from .. import fitters_ZPK
from .. import annotate
from .. import plots
from .. import TFmath
from ..utilities import args

from . import order_reduce
from . import disc_sequence
from . import disc_sequence_mag
from . import fit_aid
from . import algorithms


def data2filter_setup(
    argB=args.UNSPEC,
    data=args.UNSPEC,
    F_Hz=args.UNSPEC,
    SNR=args.UNSPEC,
    SNR_initial=args.UNSPEC,
    SNR_cutoff=args.UNSPEC,
    F_nyquist_Hz=args.UNSPEC,
    ZPK=args.UNSPEC,
    ZPK_overlay=args.UNSPEC,
    order_initial=args.UNSPEC,
    verbosity_limit=args.UNSPEC,
    delay_s=args.UNSPEC,
    aid=args.UNSPEC,
    dropnans=args.UNSPEC,
    select=args.UNSPEC,
    complete=args.UNSPEC,
    hints=None,
    doc_db=None,
):
    data = args.argscan(locals(), argB, args.REQ, arg="data")
    F_Hz = args.argscan(locals(), argB, args.REQ, arg="F_Hz")
    SNR = args.argscan(locals(), argB, 1, arg="SNR")
    SNR_initial = args.argscan(locals(), argB, None, arg="SNR_initial")
    SNR_cutoff = args.argscan(locals(), argB, 100, arg="SNR_cutoff")
    F_nyquist_Hz = args.argscan(locals(), argB, 0, arg="F_nyquist_Hz")
    ZPK = args.argscan(locals(), argB, None, arg="ZPK")
    ZPK_overlay = args.argscan(locals(), argB, None, arg="ZPK_overlay")
    order_initial = args.argscan(locals(), argB, None, arg="order_initial")
    verbosity_limit = args.argscan(locals(), argB, 5, arg="verbosity_limit")
    delay_s = args.argscan(locals(), argB, 0, arg="delay_s")
    aid = args.argscan(locals(), argB, None, arg="aid")
    dropnans = args.argscan(locals(), argB, True, arg="dropnans")
    select = args.argscan(locals(), argB, None, arg="select")
    complete = args.argscan(locals(), argB, True, arg="complete")

    if delay_s == "fit":
        delay_s = None

    if F_nyquist_Hz == 0:
        F_nyquist_Hz = None

    if aid is None:
        # TODO make annotate_setup
        doc_db = annotate.create(
            doc_db, name="fit_sequence version 1", verbosity_limit=verbosity_limit
        )
        aid = fit_aid.FitAid(
            doc_db=doc_db,
            hints=hints,
        )
    else:
        doc_db = aid.doc_db

    SNR = np.minimum(SNR, SNR_cutoff)

    if SNR_initial is None:
        if SNR is None:
            SNR_initial = 1
        else:
            SNR_initial = SNR
    else:
        SNR_initial = np.minimum(SNR_initial, SNR_cutoff)

    if select is not None:
        (F_Hz, data, SNR, SNR_initial) = np.broadcast_arrays(
            F_Hz,
            data,
            SNR,
            SNR_initial,
        )
        F_Hz = F_Hz[select]
        data = data[select]
        SNR = SNR[select]
        SNR_initial = SNR_initial[select]

    if dropnans:
        nan_select = np.isfinite(data)
        if not np.all(nan_select):
            (F_Hz, data, SNR, SNR_initial,) = np.broadcast_arrays(
                F_Hz,
                data,
                SNR,
                SNR_initial,
            )
            F_Hz = F_Hz[nan_select]
            data = data[nan_select]
            SNR = SNR[nan_select]
            SNR_initial = SNR_initial[nan_select]
    return wavestate.bunch.Bunch(locals())


def data2filter(argB=args.UNSPEC, mag_only=args.UNSPEC, alt_res=args.UNSPEC, **kwargs):
    mag_only = args.argscan(locals(), argB, False, arg="mag_only")
    alt_res = args.argscan(locals(), argB, False, arg="alt_res")
    setup = data2filter_setup(argB, **kwargs)

    aid = setup.aid

    aid.hint_set("resavg_RthreshOrdDn", 1.10, default=True)
    aid.hint_set("resavg_RthreshOrdUp", 1.00, default=True)
    aid.hint_set("resavg_RthreshOrdC", 1.00, default=True)
    aid.hint_set("overlay_native", False, default=True)

    if mag_only:
        aid.hint_set("residuals_type", "log", default=True)
        aid.hint_set("residuals_type_alt", "log", default=True)
    elif not alt_res:
        aid.hint_set("residuals_type", "log", default=True)
        aid.hint_set("residuals_type_alt", "dualA", default=True)
    else:
        aid.hint_set("residuals_type", "dualA", default=True)
        aid.hint_set("residuals_type_alt", "log", default=True)

    aid.annotate(
        name=None,
        method="v1.fit_sequence",
        about=(
            """
            Version 1 smart fitter in IIRrational library.
            Uses SVD method with high order over-fitting,
            then switches to nonlinear fits with heuristics
            to remove poles and zeros down to a reasonable system order.
            """
        ),
        verbosity=3,
    )
    MRFkwargs = dict(
        residuals_type=aid.hint("residuals_type", default="log"),
    )
    if mag_only:
        MRFkwargs["residuals_log_im_scale"] = 0
    else:
        pass

    overlay_native = aid.hint("overlay_native", default=False)

    if setup.ZPK is None:
        if not mag_only:
            ratdisc = disc_sequence.rational_disc_fit
        else:
            ratdisc = disc_sequence_mag.rational_disc_fit_mag

        if not overlay_native and setup.ZPK_overlay is not None:
            ZPK_overlay_hold = setup.ZPK_overlay
            ZPK_overlay = ((), (), 1)
            data_hold = setup.data
            data = data_hold / TFmath.TF_ZPK(
                F_Hz=setup.F_Hz,
                ZPK=ZPK_overlay_hold,
                F_nyquist_Hz=setup.F_nyquist_Hz,
            )
        else:
            data = setup.data
            ZPK_overlay = setup.ZPK_overlay

        with aid.annotate_into("initial"):
            rat_fitter = ratdisc(
                F_Hz=setup.F_Hz,
                data=data,
                SNR=setup.SNR_initial,
                nyquist_final_Hz=setup.F_nyquist_Hz,
                order=setup.order_initial,
                ZPK_overlay=ZPK_overlay,
                aid=aid,
                complete=setup.complete,
            )
        # secret output to log the rational fit output

        if not overlay_native and ZPK_overlay is not None:
            Zov, Pov, Kov = ZPK_overlay_hold
            rat_fitter.zeros_overlay = Zov
            rat_fitter.poles_overlay = Pov
            rat_fitter.gain = rat_fitter.gain * Kov
            rat_fitter.data = data_hold

        fitter = fitters_ZPK.ZPKrep2MRF(
            rat_fitter.ZPKrep,
            delay_hide=True,
            SNR=setup.SNR,
            coding_map=fitters_ZPK.coding_maps.nlFBW_safe,
            **MRFkwargs
        )

        validate = aid.hint("rational_fitter_validate", default=None)
        if validate is not None:
            validate(rat_fitter, fitter)

        # TODO, add documentation
        aid.fitter_update(fitter)
    else:
        # TODO, just use ZPKrep constructor
        aid.fitter_update(
            fitters_ZPK.ZPKdata2MRF(
                F_Hz=setup.F_Hz,
                data=setup.data,
                SNR=setup.SNR,
                ZPK=setup.ZPK,
                ZPK_overlay=setup.ZPK_overlay,
                F_nyquist_Hz=setup.F_nyquist_Hz,
                delay_hide=True,
                coding_map=fitters_ZPK.coding_maps.nlFBW_safe,
                **MRFkwargs
            )
        )

    aid.log(
        "Initial Order: (Z= {0}, P= {1}, Z-P= {2})".format(
            len(aid.fitter.zeros),
            len(aid.fitter.poles),
            len(aid.fitter.zeros) - len(aid.fitter.poles),
        )
    )

    duals = order_reduce.match_pairs(aid=aid)
    # do a pre-reduce for speed with very conservative cutoff
    ORb = order_reduce.order_reduce(
        aid=aid,
        duals=duals,
        Q_rank_cutoff=0.6,
    )
    aid.fitter_update(ORb.fitter_out)

    opt_r = aid.fitter.optimize(
        residuals_type=aid.hint("residuals_type_alt", default="log"),
        aid=aid,
    )
    aid.fitter_update()

    opt_r = aid.fitter.optimize(aid=aid)
    aid.fitter_update()

    aid.annotate(
        name="nonlinear pre-reduce and optimize",
        fitter=aid.fitter,
        plotter=plots.plot_fitter_flag,
        method="MultiReprFilterZ.optimize",
        about=(
            """
            initial conservative order reduction
            (for speed), followed by a nonlinear optimization.
            """
        ),
        verbosity=3,
    )

    ##TODO broken
    # if aid.hint('order_increase', default = True):
    #    improved = algorithms.insert_triplets(aid = aid)
    #    if improved:
    #        improved = algorithms.insert_triplets(aid = aid)

    opt_r = aid.fitter.optimize(aid=aid)
    aid.fitter_update()

    algorithms.set_min_BW(aid=aid)

    opt_r = aid.fitter.optimize(aid=aid)
    aid.fitter_update()

    aid.annotate(
        name="optimize nonlinear after bandwidth limiting",
        fitter=aid.fitter,
        plotter=plots.plot_fitter_flag,
        method="MultiReprFilterZ.optimize",
        about=(
            """
            Root bandwidths limited in the nonlinear representation
            to half of the local average the frequency spacing.
            Nonlinear optimization then applied.
            """
        ),
        verbosity=7,
    )

    with aid.annotate_into("order reduce 1"):
        order_reduce.order_reduce_aggressive(aid=aid)

    if aid.fitter.F_nyquist_Hz is not None:
        with aid.annotate_into("reduce_negroots"):
            order_reduce.order_reduce_negroots(aid=aid)
    aid.fitter_update()

    with aid.annotate_into("remove_weakroots"):
        improved = True
        while improved:
            improved = order_reduce.order_reduce_weakroots(
                aid=aid,
                delay_s=setup.delay_s,
            )
            aid.fitter_update()

    if aid.fitter.F_nyquist_Hz is not None:
        with aid.annotate_into("reduce_negroots"):
            order_reduce.order_reduce_negroots(aid=aid)

    opt_r = aid.fitter.optimize(
        aid=aid,
        residuals_type=aid.hint("residuals_type_alt", default="log"),
    )
    aid.fitter_update()

    opt_r = aid.fitter.optimize(
        aid=aid,
    )
    aid.fitter_update()

    if aid.fitter.F_nyquist_Hz is not None:
        with aid.annotate_into("reduce_negroots"):
            order_reduce.order_reduce_negroots(aid=aid)

    if setup.delay_s is None:
        aid.fitter.optimize_delay(aid=aid)
        aid.fitter_update()

    if aid.hint("flip_mindelay", default=True):
        with aid.annotate_into("flip_mindelay"):
            algorithms.flip_mindelay_opt(aid=aid)

    aid.fitter.optimize(aid=aid)
    aid.fitter_update()

    if setup.delay_s is None:
        aid.fitter.optimize_delay(aid=aid)
        aid.fitter_update()

    # final check that the sign is right..
    algorithms.sign_check_flip(aid.fitter)
    aid.fitter.optimize(aid=aid)
    aid.fitter_update()

    if setup.delay_s is None:
        aid.log("DELAY: ", aid.fitter.gain_delay_coding.p_delay)

    aid.log("FINAL RESIDUALS", aid.fitter.residuals_average)

    validate = aid.hint("fitter_final_validate", default=None)
    if validate is not None:
        validate(rat_fitter)
    # TODO, check sign and flip

    aid.annotate(
        name="final",
        fitter=aid.fitter,
        plotter=plots.plot_fitter_flag,
        method="MultiReprFilter.optimize",
        about=(
            """
            Optimizes using nonlinear parameterizations,
            """
        ),
        verbosity=0,
    )

    return aid
