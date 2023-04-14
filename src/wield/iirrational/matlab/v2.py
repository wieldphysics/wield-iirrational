#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from wield.iirrational.v2 import (
    data2filter,
    ResultsAid,
)

# from wield.iirrational.RDF import RationalDiscFilter
# from wield.iirrational.RDF2MRF import RDF2MRF
# from wield.iirrational.MRF import MultiReprFilterZ
from wield.iirrational.testing import IIRrational_data

# monkeypatch in the msurrogate annotations
ResultsAid._msurrogate_MT = False
data2filter._msurrogate_MT = True

__all__ = [
    "data2filter",
    "IIRrational_data",
]
