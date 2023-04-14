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
import collections
from wield.bunch import Bunch

from .base import (
    ArgumentError,
    mapcheck_bool,
    mapcheck_float,
    mapcheck_positive_float_orNone,
    cplx_iIjJ_list,
)

from ... import representations


def check_F_Hz(aid, hname, val):
    val = np.asarray(val)
    if np.any(val < 0):
        raise ArgumentError("Only positive frequencies allowed in 'F_Hz'")


def list_flatten(l):
    l_flat = []

    def flatten(l):
        if isinstance(l, (list, tuple)):
            for subl in l:
                flatten(subl)
        elif isinstance(l, np.ndarray):
            l_flat.extend(l)
        else:
            l_flat.append(l)

    flatten(l)
    return l_flat


def normalize_roots(val):
    if isinstance(val, str):
        if val.lower() == "none":
            val = None
        else:
            val = cplx_iIjJ_list(val)
    elif isinstance(val, collections.Mapping):
        r = val["real"]
        i = val["imag"]
        val = np.asarray(r) + 1j * np.asarray(i)
    if val is None:
        return None
    # TODO, treat the single number case better
    val = list_flatten([val])
    return val


def mapcheck_roots(aid, hname, val):
    try:
        val = normalize_roots(val)
        complete = aid.hint("pair_completion")
        return representations.asMRRB(val, complete=complete)
    except ValueError:
        raise ArgumentError(
            ("Argument {} must be a list of complex numbers").format(hname, val)
        )


def mapcheck_ZPK(aid, hname, val):
    if val is None:
        return None
    complete = aid.hint("pair_completion")
    return representations.asZPKTF(val, complete=complete)


def mapcheck_mode(aid, hname, val):
    modes = [
        "full",
        "fullAAA",
        "full2x",
        "fit",
        "AAA",
        "onlyAAA",
        "onlyAAAreduce",
        "reduce",
        "rational",
        "rational2x",
        "dumpargs",
        "chebydebug",
        "chebydebug+",
        "discdebug",
        "copy",
        "gain",
    ]
    if val not in modes:
        raise ArgumentError("Argument {} must be one of {}".format(hname, modes))
    return val


kw_hints_pre = Bunch(
    pair_completion=dict(
        mapcheck=mapcheck_bool,
        default=True,
    )
)


kw_hints_data = Bunch(
    data=dict(
        APgroup="fitting",
        # TODO, mapcheck
        aliases=["xfer"],
        aliases_bad=["TF", "d", "transfer"],
        about="""
        Data of the transfer function to fit. Should be complex. Can only fit
        transfer functions on real data (conjugate negative frequencies).
        """,
    ),
    F_Hz=dict(
        APgroup="fitting",
        # TODO, mapcheck
        aliases=["f_Hz", "frequency_Hz"],
        aliases_bad=["freq", "F", "f", "frequency"],
        about="""
        Frequencies for the transfer function data. Should only be positive.
        """,
        checkmap=check_F_Hz,
    ),
    SNR=dict(
        APgroup="fitting",
        # TODO, mapcheck
        aliases=["W", "snr"],
        aliases_bad=["weight", "weights"],
        default=None,
        about="""
        Statistical weights used for fitting and decision procedures.
        Values not trusted and rescaled adjusted in default algorithms. Use
        the 'emphasis' parameter to adjust weighting between statistical
        decisions and accuracy/implementation requirements for the fit.
        """,
    ),
    coherence=dict(
        APgroup="fitting",
        # TODO, mapcheck
        aliases=["coherence", "coh", "COH"],
        default=None,
        about="""
        Coherence of transfer functions used. This is converted into SNR.
        """,
    ),
    SNR_phase_relative=dict(
        APgroup="fitting",
        # TODO, mapcheck
        aliases=["Wphase_relative", "SNR_phase_rel"],
        default=None,
        about="""
        Statistical weights specifically for the phase part of the fit. This
        only affects the ZPK fitting (second half) of "full" mode, and the
        "fit", "reduce" and "gain" modes. The default of "None" specifies to
        use the SNR for both types of fit. This is a relative scaling to adjust
        the SNR by for the phase aspect of the fit.
        """,
    ),
    emphasis=dict(
        # TODO, mapcheck
        APgroup="fitting",
        aliases=["Wemph"],
        default=None,
        about="""
        Emphasis weighting. Affects fits and some decision procedures differently
        than the purely statistical data. Not rescaled or adjusted by algorithms.
        Use values with small dynamic range as the reduction in the effective
        number of data points scales with the dynamic range of the emphasis squared.
        """,
    ),
    select=dict(
        APignore=True,
        # TODO, mapcheck
        default=None,
        about="""
        This may be a slice object or boolean array to remove points from: data,
        F_Hz, SNR, and emphasis. May be anything that can index numpy arrays.
        """,
    ),
)

kw_hints_fit = Bunch(
    help=dict(
        APignore=True,
        aliases=["h"],
        mapcheck=mapcheck_bool,
        default=False,
        about="""
        Print a help screen
        """,
    ),
    # TODO, added new modes to docs
    mode=dict(
        APgroup="fitting",
        APpriority=3,
        mapcheck=mapcheck_mode,
        default="full",
        about="""
        Fitting mode to use, to change the automation level. Must be one of
         - "full": to use rational fitting to get initial parameters, then alternate
                   optimization and order reduction
         - "rational": to use rational fitting to get initial parameters and then
                       stop after the simplest order reduction. Useful to use with
                       other fitters. Does not perform delay fitting.
         - "fit": To take an initial ZPK guess, and only fit/optimize it. Good
                  for refining previous fits.
         - "reduce": To take an initial ZPK guess, and then alternate fitting and
                     order reduction.
         - "dumpargs": Return dictionary of arguments to see the full run settings.
        """,
    ),
    chain=dict(
        APgroup="fitting",
        APignore=True,
        APpriority=3,
        # mapcheck = mapcheck_chain,
        default=None,
        about="""
        Chain argument to run sequential fits. Takes a fitter results object.
        This will change the default fitting from mode='full' to
        mode='fit_reduce', but will keep mode='fit'. This jumpstarts all
        defaults to use those supplied for the previous fit. In addition, the
        ZPKs of the previous fit are used. Any additional arguments will override
        the defaults. This can be used to modify the fitting SNR/weights, or
        possibly to activate H_infinity weighting.
        """,
    ),
    F_nyquist_Hz=dict(
        APshort="-N",
        APignore=False,
        # TODO, mapcheck
        APpriority=5,
        mapcheck=mapcheck_positive_float_orNone,
        default=None,
        about="""
        This selects the Nyquist frequency for Z-domain fitting. If None or not
        specified, the fits are done in the "Sf" domain - the s domain scaled for
        units of frequency.
        """,
    ),
    ZPK=dict(
        APignore=True,
        mapcheck=mapcheck_ZPK,
        require_hints=["pair_completion"],
        aliases=["zpk"],
        default=None,
    ),
    zeros=dict(
        APshort="-z",
        APnargs="*",
        APgroup="fitting",
        APtype=cplx_iIjJ_list,
        mapcheck=mapcheck_roots,
        normalize=normalize_roots,
        require_hints=["pair_completion"],
        aliases=["Z", "z"],
        default=None,
        about="Initial zeros used in the fit.",
    ),
    poles=dict(
        APshort="-p",
        APgroup="fitting",
        APnargs="*",
        APtype=cplx_iIjJ_list,
        mapcheck=mapcheck_roots,
        normalize=normalize_roots,
        require_hints=["pair_completion"],
        aliases=["P", "p"],
        default=None,
        about="Initial poles used in the fit.",
    ),
    gain=dict(
        APshort="-k",
        APgroup="fitting",
        APtype=float,
        mapcheck=mapcheck_float,
        aliases=["K", "k"],
        default=None,
        about="Initial gain for the fit",
    ),
    ZPK_overlay=dict(
        APignore=True,
        mapcheck=mapcheck_ZPK,
        require_hints=["pair_completion"],
        aliases=["zpk_overlay"],
        default=None,
    ),
    zeros_overlay=dict(
        APshort="-Z",
        APgroup="fitting",
        APnargs="*",
        APtype=cplx_iIjJ_list,
        mapcheck=mapcheck_roots,
        normalize=normalize_roots,
        require_hints=["pair_completion"],
        aliases=["Zov"],
        default=None,
        about="zeros guaranteed to be part of the fit.",
    ),
    poles_overlay=dict(
        APshort="-P",
        APgroup="fitting",
        APnargs="*",
        APtype=cplx_iIjJ_list,
        mapcheck=mapcheck_roots,
        normalize=normalize_roots,
        require_hints=["pair_completion"],
        aliases=["Pov"],
        default=None,
        about="poles guaranteed to be part of the fit.",
    ),
)

kw_hints = dict()
kw_hints.update(kw_hints_fit)
kw_hints.update(kw_hints_data)
