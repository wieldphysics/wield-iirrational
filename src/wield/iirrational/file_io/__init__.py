#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from .any_io import (
    load_any,
    write_any,
)

from .csv_io import (
    load_csv,
)

from .types import (
    type2type,
    ext2type,
    type2features,
    determine_type,
)

from .utilities import (
    subkey_search,
)


def save(fname, d):
    typeB = determine_type(fname)
    write_any(
        fname=typeB.fname,
        ftype=typeB.ftype,
        fdict=d,
    )


def load(fname, ftype=None):
    typeB = determine_type(fname)
    if ftype is None:
        ftype = typeB.ftype

    return load_any(
        fname=typeB.fname,
        ftype=ftype,
    )
