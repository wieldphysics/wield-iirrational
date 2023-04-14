#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from wield.iirrational.testing.utilities import (
    sign_validate_and_plot_hint,
    sign_validate_and_digest_hint,
)

from ..testing.plots import plot_on_assert


def validate_plot_log(fname, request):
    hdict = sign_validate_and_plot_hint(__file__, request)
    hdict.update(
        dict(
            log_level_debug=10,
            log_level=10,
            log_level_warn=10,
        )
    )
    return hdict
