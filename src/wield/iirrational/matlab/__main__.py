#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import os
from msurrogate import SurrogateApp

from . import v1
from . import testing

from ..plots import IIRPlots

plots = IIRPlots(
    show_immediate=True,
    show_saved=False,
    force_save=False,
)

if __name__ == "__main__":
    app = SurrogateApp()
    # app.option('serializer-preferred', 'dill')
    app.option("serializer-preferred", "pickle")
    metaD = app.meta_daemon_setup()
    metaD.register(v1, "v1")
    metaD.register(v1, "v2")
    metaD.register(testing, "testing")
    metaD.register(plots, "plots")
    app.run_loop()
