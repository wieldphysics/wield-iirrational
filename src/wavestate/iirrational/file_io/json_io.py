#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import json
import sys


def load_json(fname):
    with open(fname) as F:
        fdict = json.load(F)
    return fdict


def write_json(fname, fdict):
    if sys.version_info < (3, 4):
        with open(fname, "w") as F:
            json.dump(fdict, F, indent=4, ensure_ascii=False)
    else:
        with open(fname, "w", encoding="utf8") as F:
            json.dump(fdict, F, indent=4, ensure_ascii=False)
    return
