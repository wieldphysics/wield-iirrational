#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import h5py
from wavestate.bunch.hdf_deep_bunch import HDFDeepBunch


def load_hdf5(fname):
    # with h5py.File(fname) as h5F:
    h5F = h5py.File(fname, "r")
    fdict = HDFDeepBunch(h5F, writeable=False)
    return fdict


def write_hdf5(fname, fdict):
    with h5py.File(fname, "w") as h5F:
        hdf = HDFDeepBunch(h5F, writeable=True)
        hdf.update_recursive(fdict)
    return
