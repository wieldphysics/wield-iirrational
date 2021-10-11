# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import declarative

from .base import (
    ArgumentError,
    mapcheck_bool,
    mapcheck_positive_int,
    mapcheck_positive_int_orNone,
)


def mc_log_level(aid, aname, val):
    val = int(val)
    if val < 0 or val > 10:
        raise ArgumentError("log_level values must be between 0-10, 0 being least verbose, 10 most.")
    return val


def hint_log_level(**kwargs):
    kw = dict(
        mapcheck = mc_log_level,
    )
    kw.update(**kwargs)
    return kw


kw_hints = declarative.Bunch(
    logging_module_use = dict(
        APignore = True,
        mapcheck = mapcheck_bool,
        about = "Use the python logging module instead of print",
        aliases = ['use_logging_module'],
    ),
    log_level = hint_log_level(
        APshort = '-l',
        APpriority = 20,
        APgroup = 'logging',
        mapcheck = mapcheck_positive_int,
        default = 5,
        about = "Log level default for all logging types",
        aliases = ['verbosity'],
    ),
    log_level_alert = hint_log_level(
        APgroup = 'logging',
        APpriority = 25,
        mapcheck = mapcheck_positive_int_orNone,
        about = """
        Data on useful statistical tests, particularly those which cause input
        data/SNR/emphasis to be reinterpreted. Typically logs at level 2-4.
        """,
        aliases = ['verbosity_alert'],
    ),
    log_level_info = hint_log_level(
        APgroup = 'logging',
        APpriority = 25,
        mapcheck = mapcheck_positive_int_orNone,
        about = """
        Logging on miscellaneous details.
        """,
        aliases = ['verbosity_info'],
    ),
    log_level_debug = hint_log_level(
        APgroup = 'logging',
        APpriority = 25,
        mapcheck = mapcheck_positive_int_orNone,
        default = 0,
        about = """
        Debugging reports used for development.
        """,
        aliases = ['verbosity_debug'],
    ),
    log_level_warn = hint_log_level(
        APgroup = 'logging',
        APpriority = 25,
        mapcheck = mapcheck_positive_int_orNone,
        default = 5,
        about = """
        Warnings about failed tests, or data operating in a regime unexpected to
        work or that validation should be made.
        """,
        aliases = ['verbosity_warn'],
    ),
    log_level_rationale = hint_log_level(
        APgroup = 'logging',
        APpriority = 25,
        mapcheck = mapcheck_positive_int_orNone,
        default = 0,
        about = """
        Detailed explanations of tests and algorithms. Extremely verbose and
        intended for first users and digest reports.
        """,
        aliases = ['verbosity_rationale'],
    ),
)

