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


def load_csv(fname, parse_str, parse_map):
    """
    Can raise KeyError if parse_str does not map into parse_map
    """
    delimiter = parse_str[0]
    if delimiter in ['-']:
        delimiter = None
    parse_str = parse_str[1:]

    if delimiter is None:
        subkeyspl = parse_str.split()
    else:
        subkeyspl = parse_str.split(delimiter)

    key_gen = []
    for sk in subkeyspl:
        key_gen.append(parse_map[sk])

    farr = np.genfromtxt(
        fname,
        delimiter = delimiter,
        filling_values = float('NaN'),
    )

    fdict = dict()
    for idx_col in range(farr.shape[0]):
        fdict[str(idx_col)] = farr[idx_col]

        if idx_col < len(subkeyspl):
            sk = key_gen[idx_col]
            if sk is not None:
                fdict[sk] = farr[idx_col]

    return fdict


