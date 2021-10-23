#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

from wavestate import declarative

from ... import fitters_ZPK
from .base import (
    ArgumentError,
    mapcheck_bool,
    mapcheck_positive_float,
    mapcheck_nonnegative_float,
    mapcheck_nonnegative_float_orNone,
    mapcheck_positive_int_orNone,
    mapcheck_nonnegative_int_orNone,
    mapcheck_positive_float_orNone,
    mapcheck_positive_int,
)


def mapcheck_downsample_type(aid, aname, val):
    val = str(val).lower()
    allowed = ["linear", "lin", "log", "linlog", "loglin"]
    if val not in allowed:
        raise ArgumentError("Argument '{}' must be one of {}".format(aname, allowed))
    return val


def mapcheck_residuals_type(aid, aname, val):
    val = str(val)
    allowed = ["log", "dualA", "dualB", "zeros", "poles"]
    if val not in allowed:
        raise ArgumentError("Argument '{}' must be one of {}".format(aname, allowed))
    return val


def default_residuals_type(aid, aname):
    if not aid.hint("alternate_residuals"):
        return "log"
    else:
        return "dualA"


def default_residuals_type_alt(aid, aname):
    if not aid.hint("alternate_residuals"):
        return "dualA"
    else:
        return "log"


def mapcheck_boostlist(aid, aname, val):
    if val is None:
        return None
    try:
        if isinstance(val, str):
            val_split = val.split(",")
            val = []
            for v in val_split:
                a, b = v.split("-")
                a = float(a)
                b = float(b)

            val.append((float(a), float(b)))
    except Exception:
        raise ArgumentError("{} must be a range tuple or list of range tuples")

    if not isinstance(val, (list, tuple)):
        raise ArgumentError("{} must be a range tuple or list of range tuples")

    for v in val:
        if not isinstance(v, (list, tuple)):
            break
        if len(v) != 2:
            raise ArgumentError("{} must be a range tuple or list of range tuples")
    else:
        return val

    if len(val) != 2:
        raise ArgumentError("{} must be a range tuple or list of range tuples")

    return (val,)


kw_hints = wavestate.bunch.Bunch(
    order_max=dict(
        APgroup="order",
        mapcheck=mapcheck_positive_int,
        default=100,
        about="""
        Maximium order that the rational fitter order estimation may use.
        Increasing this value from the default may impact numerical stability.
        """,
    ),
    order_min=dict(
        APgroup="order",
        mapcheck=mapcheck_positive_int,
        default=20,
        about="""
        Minimum order to use during order estimation. Smaller values will speed
        up fits, but may fail to fix complex data.
        """,
    ),
    order_initial=dict(
        APgroup="order",
        mapcheck=mapcheck_positive_int_orNone,
        default=None,
        aliases=[
            "order",
            "rational_fit_order",
            "rational_cheby_fit_order",
        ],
        about="""
        Order to use in rational fitting. If not specified, the order is
        increased until the residuals plateau.
        """,
    ),
    order_first=dict(
        APgroup="order",
        mapcheck=mapcheck_positive_int_orNone,
        default=8,
        aliases=[
            "order1",
        ],
        about="""
        Order to use in rational fitting for the first round of full2x. Can be lower to be fast.
        """,
    ),
    # emphasis_order_max = dict(
    #    APgroup = 'emphasis',
    #    mapcheck = mapcheck_positive_int_orNone,
    #    default = None,
    #    about = """
    #    Maximium order that the rational fitter order estimation may use during
    #    the emphasis application stage.
    #    """,
    # ),
    # emphasis_order_min = dict(
    #    APgroup = 'emphasis',
    #    #TODO, check minmax with emphasis_order_max
    #    mapcheck = mapcheck_positive_int_orNone,
    #    default = None,
    #    about = """
    #    Minimum order to use during order estimation during the emphasis application
    #    stage.
    #    """,
    # ),
    # emphasis_order_initial = dict(
    #    APgroup = 'emphasis',
    #    mapcheck = mapcheck_positive_int_orNone,
    #    default = None,
    #    aliases = [
    #        'emphasis_order',
    #        'emphasis_rational_fit_order',
    #        'emphasis_rational_cheby_fit_order',
    #    ],
    #    about = """
    #    Order to use in rational fitting during the emphasis application stage.
    #    """,
    # ),
    F_boost_Hz=dict(
        APgroup="SNR adjustments",
        mapcheck=mapcheck_boostlist,
        default=None,
        about="""
        start-end frequency pairs (e.g. 1.2-1.6) to double the SNR end to add emphasis during fitting. Multiple can be specified, comma separated.
        """,
    ),
    downsample=dict(
        APgroup="SNR adjustments",
        mapcheck=mapcheck_positive_int_orNone,
        default=None,
        about="""
        Downsample the data to this number of points. Uses SNR weighing to average during downsampling.
        """,
    ),
    downsample_type=dict(
        APgroup="SNR adjustments",
        default="log",
        about="""
        Use this type of spacing for the downsampled points: linear, log, or loglinear.
        """,
    ),
    F_max_Hz=dict(
        APgroup="SNR adjustments",
        mapcheck=mapcheck_positive_float_orNone,
        default=None,
        about="""
        Maximum frequency to use, cuts off data above this frequency.
        """,
    ),
    F_min_Hz=dict(
        APgroup="SNR adjustments",
        mapcheck=mapcheck_positive_float_orNone,
        default=None,
        about="""
        Minimum frequency to use, cuts off data below this frequency.
        """,
    ),
    SNR_max=dict(
        APgroup="SNR Adjustments",
        mapcheck=mapcheck_positive_float_orNone,
        default=None,
        aliases=[
            "W_max",
        ],
        about="""
        Maximum SNR. Hard cutoff for SNR above this level.
        """,
    ),
    SNR_min=dict(
        APgroup="SNR Adjustments",
        mapcheck=mapcheck_nonnegative_float_orNone,
        default=None,
        aliases=[
            "W_min",
        ],
        about="""
        Minimum SNR. Hard cutoff for SNR below this level.
        """,
    ),
    SNR_regularize_ratio=dict(
        APgroup="SNR Adjustments",
        mapcheck=mapcheck_positive_float,
        default=0.5,
        about="""
        Ratio of effective data points, determined by the dynamic range of the
        SNR weighting, to the actual number of data points. This parameter
        adjusts how the SNR is adjusted to ensure sufficient data points for a
        good fit.
        """,
    ),
    SNR_regularize_scale=dict(
        APgroup="SNR Adjustments",
        mapcheck=mapcheck_positive_float,
        default=10,
        about="""
        Similar to SNR_regularize_ratio, but instead of a direct ratio, the ratio
        is determined as 'SNR_regularize_scale' / max(SNR). For very high SNR
        fits, the SNR should be whitened significantly using the ratio adjustment,
        since the error is dominated by systematics rather than statistical noise.
        The default is 10, causing a %90 coverage ratio when the highest SNR is
        100 (1% statistical error).
        """,
    ),
    SNR_estimate_width=dict(
        APgroup="SNR Adjustments",
        mapcheck=mapcheck_nonnegative_int_orNone,
        default=10,
        about="""
        The window width of nearby points to use for averaging and median filtering
        when estimating the SNR via windowed sample variance. Default is 10.
        None or 0 indicates to not try.
        """,
    ),
    trust_SNR=dict(
        APgroup="SNR Adjustments",
        APaction="store_true",
        mapcheck=mapcheck_bool,
        default=False,
        about="""
        Overall parameter determining if SNR/W statistical weighting input should
        be trusted. If False (the default) then additional statistical tests are
        run to adjust the SNR for better fits. If the user can fully trust their
        SNR estimators and fits to be unbiased, then the user likely does not
        need this tool.
        """,
    ),
    inverse=dict(
        APgroup="fit",
        APaction="store_true",
        mapcheck=mapcheck_bool,
        default=False,
        about="""
        take the reciprocal/inverse of the data (1/xfer) and fit that. Good for fitting compensation filters.
        """,
    ),
    never_unstable_poles=dict(
        APgroup="fit",
        APaction="store_true",
        mapcheck=mapcheck_bool,
        default=False,
        about="""
        During the phase patching tests, unstable roots will be added if detected
        with statistical significance. This prevents the addition of unstable poles.
        """,
    ),
    never_unstable_zeros=dict(
        APgroup="fit",
        APaction="store_true",
        mapcheck=mapcheck_bool,
        default=False,
        about="""
        During the phase patching tests, unstable roots will be added if detected
        with statistical significance. This prevents the addition of unstable zeros.
        """,
    ),
    suggest=dict(
        APgroup="operating mode",
        APaction="store_true",
        mapcheck=mapcheck_bool,
        default=False,
        about="""
        How to use the rational fitter when provided an initial ZPK argument.
        the default of False causes the rational to overlay the ZPK suggestion,
        forcing the roots in the suggestion to be used during the optimization
        stage. If set to True, the ZPK will suggest initial poles to use during
        the fit, which can accellerate convergence or cause it to start at a
        lower initial order, speeding up future fitting operations.
        """,
    ),
    coding_map=dict(
        APignore=True,
        # TODO, add a mapcheck
        default=fitters_ZPK.codings_s.coding_maps.SOSnl,
        advanced=True,
        about="""
        The coding map for the nonlinear fitter. This determines the nonlinear
        parameterization best suited for the NLS fine-grained optimization. These
        parameterizations can also influence the bounds and tradoffs of roots
        placement.
        """,
    ),
    resavg_RthreshOrdDn=dict(
        APgroup="advanced",
        default=1.10,
        mapcheck=mapcheck_positive_float,
        about="""
        The threshold relative change in the average residuals to accept a fit
        of lower order. Only used for the "baseline" fit determination before
        delay is activated and not used during the total order reduction.
        """,
    ),
    resavg_RthreshOrdUp=dict(
        APgroup="advanced",
        default=0.95,
        mapcheck=mapcheck_positive_float,
        about="""
        The threshold relative change in the average residuals to accept a fit
        of higher order. Only used for the "baseline" fit determination before
        delay is activated and not used during the total order reduction.
        """,
    ),
    resavg_RthreshOrdC=dict(
        APgroup="advanced",
        default=1.00,
        mapcheck=mapcheck_positive_float,
        about="""
        The threshold relative change in the average residuals to accept a fit
        of equal order. Only used for the "baseline" fit determination before
        delay is activated and not used during the total order reduction.
        """,
    ),
    distance_limit_scale=dict(
        APgroup="tunings",
        default=1,
        mapcheck=mapcheck_positive_float,
        about="""
        Scaling for how to limit root bandwidth for roots in-between data points.
        Increase to force roots to lower Q, decrease to allow higher Qs.
        """,
    ),
    alternate_residuals=dict(
        APgroup="residuals",
        default=False,
        mapcheck=mapcheck_bool,
        about="""
        Flag to switch from using log/phase residuals to using dual
        ratio residuals as the principle fit residuals. Log/phase is typically
        better behaved.
        """,
    ),
    h_infinity=dict(
        APgroup="residuals",
        # TODO, may a checker for this
        # mapcheck = mapcheck_nonnegative_float,
        default=0,
        about="""
        Use an partial h_infinity norm during the nonlinear fitting
        (does not apply to the linear fits). This argument takes a [0, 1] float
        for the fraction of the residuals to include. If 1, this is true h_infinity,
        if < 1, then the residuals are ranked, and the smallest fraction of them
        have h_infinity_deweight applied to them, to lower their affect of the fit.
        This argument may also be an array, in which case it reweights the
        ranked residuals (ranking array should be largest to smallest). It may
        alternatively take a function, which is given the length, and returns
        the reweighting array.
        """,
    ),
    h_infinity_deweight=dict(
        APgroup="residuals",
        mapcheck=mapcheck_nonnegative_float,
        default=0.1,
        about="""
        Factor to de-weight the residuals below the h_infinity threshold. If 0,
        they are clipped. Defaults to 0.1, to reduce the significance 10x.
        Factors >1 imply an h_0 type fit, where outliers are more strongly ignored.
        """,
    ),
    residuals_type=dict(
        APgroup="residuals",
        APchoices=["log", "dualA", "dualB", "poles", "zeros"],
        default=default_residuals_type,
        mapcheck=mapcheck_residuals_type,
        advanced=True,
        require_hints=["alternate_residuals"],
        about="""
        Standard residuals type to use for the optimizations. Must be one of the
        definitions below, where R=fit/data, the ratio of fit to data and W is
        the SNR:
            log:   W*(ln(|R|) + 1j*R.imag/|R|)
            dualA: W*(R + 1/R - 2)/2
            dualB: W*(R - 1/R)/2
            poles: W*(1/R - 1)
            zeros: W*(R - 1)
        """,
    ),
    residuals_type_alt=dict(
        APgroup="residuals",
        default=default_residuals_type_alt,
        APchoices=["log", "dualA", "dualB", "poles", "zeros"],
        mapcheck=mapcheck_residuals_type,
        advanced=True,
        require_hints=["alternate_residuals"],
        about="""
        Standard alternate residuals type to use during optimization annealing.
        See residuals_type for options.
        """,
    ),
    prune_Qrank=dict(
        APgroup="tuning",
        default=None,
        mapcheck=mapcheck_positive_float_orNone,
        aliases=[
            "prune",
        ],
        about="""
        Pre-prune the filter. Useful for order reduction of StateSpace models.
        The Q-ranking is like a weighted fminreal order reduction. It ranks
        pole-zero by their distance divided by bandwidth. Values smaller than
        .1-.5 remove only pairs that have small adjustments to a transfer function
        and can generally be removed without modifying it.
        """,
    ),
    root_bandwidth_Hz_max=dict(
        APgroup="tunings",
        default=None,
        mapcheck=mapcheck_positive_float,
        aliases=[
            "max_BW_Hz",
            "BW_max_Hz",
        ],
        aliases_bad=[
            "root_bandwidth_max",
            "BW_max",
            "max_BW",
        ],
        about="""
        Maximum bandwidth of any pole or zero in the fit. This affects the
        asymptotic rolloff when relative_degree is constrained.
        """,
    ),
    root_F_Hz_max=dict(
        APgroup="tunings",
        default=None,
        mapcheck=mapcheck_positive_float,
        aliases=[
            "max_F_Hz",
        ],
        aliases_bad=[
            "root_F_max",
            "F_max",
            "max_F",
        ],
        about="""
        Maximum frequency of root complex part.
        """,
    ),
    greedy_order=dict(
        APgroup="speed",
        default=30,
        mapcheck=mapcheck_positive_int_orNone,
        about="""
        Does only partial optimization during the reduce step. Greatly speeds
        up order reduction, but may not be as effective. Defaults to 30, where
        it uses greedy optimization until reaching order 30, then uses combinatoric
        optimization.
        """,
    ),
    baseline_only=dict(
        APgroup="speed",
        default=False,
        APaction="store_true",
        mapcheck=mapcheck_bool,
        about="""
        Only fit down to baseline, do not do successive order reduction.
        """,
    ),
    multithreading=dict(
        APgroup="computing",
        default=None,
        mapcheck=mapcheck_positive_int_orNone,
        about="""
        Perform the order reduction trials using a thread pool of this size.
        If None (the default) or 1, do not activate the pool. Note that many python/numpy
        implementations use a BLAS library that already multithread linear algebra operations.
        """,
    ),
)
