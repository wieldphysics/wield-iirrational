#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import pickle


def load_pickle(fname):
    with open(fname) as F:
        fdict = pickle.load(F)
    return fdict


def write_pickle(fname, fdict):
    with open(fname, "wb") as F:
        fdict = pickle.dump(fdict, F)
    return fdict
