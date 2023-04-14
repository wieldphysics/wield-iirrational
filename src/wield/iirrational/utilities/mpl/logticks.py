#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from matplotlib import ticker
import numpy as np


def ticks_log_format(value, index):
    """
    get the value and returns the value as:
       integer: [0,99]
       1 digit float: [0.1, 0.99]
       n*10^m: otherwise
    To have all the number of the same size they are all returned as latex strings
    """
    exp = np.floor(np.log10(value))
    base = value / 10 ** exp
    if exp == 0 or exp == 1:
        return "${0:d}$".format(int(value))
    if exp == -1:
        return "${0:.1f}$".format(value)
    else:
        return "${0:d}\\times10^{{{1:d}}}$".format(int(base), int(exp))


def setup_log_xticks(ax):
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(ticks_log_format)
    )  # add the custom ticks
    ax.xaxis.set_minor_formatter(
        ticker.FuncFormatter(ticks_log_format)
    )  # add the custom ticks
