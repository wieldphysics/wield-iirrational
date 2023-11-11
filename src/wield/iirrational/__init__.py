#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from os import path

from .version import (
    version,
    __version__,
)

from . import v1
from . import v2
from .file_io import save, load


def matlabpath():
    # import msurrogate
    # don't join them as it only seems to work on linux
    return ":".join(
        [
            # msurrogate.matlabpath(),
            path.abspath(path.split(__file__)[0]),
        ]
    )


__all__ = [
    "version",
    "__version__",
    "v1,"
    "v2",
    "save,"
    "load",
    "matlabpath,",
    "fitters_ZPK"
]
