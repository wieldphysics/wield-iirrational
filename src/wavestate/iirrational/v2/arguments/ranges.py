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

from .base import (
    ArgumentError,
    mapcheck_positive_float_orNone,
    mapcheck_positive_float,
    mapcheck_int_orNone,
    mapcheck_positive_int_orNone,
)

def mapcheck_reldeg(aid, aname, val):
    val = mapcheck_int_orNone(aid, aname, val)
    if val is None:
        return val
    ival = int(val)
    if ival != float(val):
        raise ArgumentError("argument '{}' must be an integer".format(aname))
    if aname == 'relative_degree_min':
        rdmax = aid.hint('relative_degree_max')
        if rdmax is not None and val > rdmax:
            raise ArgumentError(
                "argument 'relative_degree_min'={} must be smaller"
                " or equal to 'relative_degree_max'={}".format(val, rdmax))
    return val

def default_reldeg_minmax(aid, aname):
    val = aid.hint('relative_degree')
    return val

def default_delay_s(aid, aname):
    val_min = aid.hint('delay_s_min')
    val_max = aid.hint('delay_s_max')
    if val_max is not None:
        if val_min < 0:
            if val_max > 0:
                val = 0
            else:
                val = val_max
        else:
            val = val_min
    else:
        val = val_min
    return val

def mapcheck_delay_s(aid, aname, val):
    val = float(val)
    val_min = aid.hint('delay_s_min')
    val_max = aid.hint('delay_s_max')
    if val < val_min:
        raise ArgumentError((
            "Argument 'delay_s'={} must be larger than 'delay_s_min'={}"
        ).format(val, val_min))

    if val_max is not None:
        if val > val_max:
            raise ArgumentError((
                "Argument 'delay_s'={} must be larger than 'delay_s_max'={}"
            ).format(val, val_max))
    elif val > val_min:
            raise ArgumentError((
                "Argument 'delay_s' must be equal to 'delay_s_min' since"
                " 'delay_s_max' is not specified. Either specify 'delay_s_max'"
                " or don't specify 'delay_s', in which case it will default to 'delay_s_min'."
            ).format(val, val_min))
    return val

def mapcheck_total_degree_min(aid, aname, val):
    val = mapcheck_positive_int_orNone(aid, aname, val)
    if val is None:
        return val
    val = int(val)
    if val < 0:
        raise ArgumentError(
            "Argument '{}'=val must not be negative".format(aname, val)
        )
    return val


kw_hints = wavestate.bunch.Bunch(
    relative_degree = dict(
        APgroup = 'order',
        APpriority = 10,
        mapcheck = mapcheck_reldeg,
        default = None,
        about = """
        Sets the initial relative degree (number zeros minus polse).
        Defaults to None, which will be the midpoint of the min and max if they
        are specified, otherwise this will default to 0, which typically will
        still fit filters to the correct degree. If constrained to be different
        than the data, the filter will enter the asymptotic regime set by the
        relative degree within 2x of the 'root_bandwidth_Hz_max' setting.
        """,
    ),
    relative_degree_max = dict(
        APgroup = 'order',
        APpriority = 11,
        mapcheck = mapcheck_reldeg,
        default = default_reldeg_minmax,
        require_hints = ['relative_degree'],
        about = """
        Maximum value for the filter relative degree (number zeros minus polse).
        Defaults to the relative degree (which may be None). If this value is
        None, the degree is unconstrained and the fit will land at some degree that fits well.
        """,
    ),
    relative_degree_min = dict(
        APgroup = 'order',
        APpriority = 12,
        mapcheck = mapcheck_reldeg,
        default = default_reldeg_minmax,
        require_hints = [
            'relative_degree',
            #needed for the ordering check
            'relative_degree_max',
        ],
        about = """
        Minimum value for the filter relative degree (number zeros minus polse).
        Defaults to the relative degree (which may be None). If this value is
        None, the degree is unconstrained and the fit will land at some degree that fits well.
        """,
    ),
    total_degree_min = dict(
        APgroup = 'order',
        APpriority = 13,
        mapcheck = mapcheck_total_degree_min,
        default = 2,
        aliases = [
            'degree_min',
        ],
        about = """
        Minimum degree to search through during the successive order reduction phase.
        Defaults to 2. If None, then successive reduction will not be performed.
        """,
    ),
    delay_s = dict(
        APgroup = 'delay',
        APpriority = 16,
        mapcheck    = mapcheck_delay_s,
        default     = default_delay_s,
        aliases     = ['delay_seconds'],
        aliases_bad = ['delay'],
        require_hints = [
            'delay_s_min',
            'delay_s_max',
        ],
        about = """
        Use this delay (in seconds) for the initial fitting up until the
        "baseline" fit determination is complete. By default this is the same
        as delay_s_min, unless delay_s_min is negative, in which case it defaults
        to the smaller of 0 seconds or delay_s_max.
        """,
    ),
    delay_s_max = dict(
        APgroup = 'delay',
        APpriority = 15,
        mapcheck = mapcheck_positive_float_orNone,
        default = None,
        aliases       = ['delay_max_s'],
        aliases_bad   = ['delay_max'],
        about = """
        The maximum delay in seconds. Defaults to None, in which case delay is
        never fit as a free parameter.
        """,
    ),
    delay_s_min = dict(
        APgroup = 'delay',
        APpriority = 17,
        mapcheck = mapcheck_positive_float,
        default = 0,
        aliases       = ['delay_min_s'],
        aliases_bad   = ['delay_min'],
        about = """
        The minimum delay in seconds. Defaults to 0. Also sets the typical
        default value of the delay used for the initial part of the fits.
        Must always be specified (cannot be None).
        """,
    ),
)

