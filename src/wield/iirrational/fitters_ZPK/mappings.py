#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from wield.bunch import Bunch
from . import codings_s
from . import codings_z

_common = dict()

coding_maps = Bunch(
    RI=Bunch(
        rep=Bunch(
            Sf=codings_s.coding_maps.RI,
            Z=codings_z.coding_maps.RI,
        ),
        **_common
    ),
    FBW=Bunch(
        rep=Bunch(
            Sf=codings_s.coding_maps.FBW,
            Z=codings_z.coding_maps.FBW,
        ),
        **_common
    ),
    nlFBW=Bunch(
        rep=Bunch(
            Sf=codings_s.coding_maps.nlFBW,
            Z=codings_z.coding_maps.nlFBW,
        ),
        **_common
    ),
    nlFBW_safe=Bunch(
        rep=Bunch(
            Sf=codings_s.coding_maps.nlFBW_safe,
            Z=codings_z.coding_maps.nlFBW_safe,
        ),
        **_common
    ),
    SOS=Bunch(
        rep=Bunch(
            Sf=codings_s.coding_maps.SOS,
            Z=codings_z.coding_maps.SOS,
        ),
        **_common
    ),
    # SOS_safe = Bunch(
    #    rep = Bunch(
    #        Sf = codings_s.coding_maps.SOS_safe,
    #        Z  = codings_z.coding_maps.SOS_safe,
    #    ),
    #    **_common
    # ),
)
