#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from .statespace import (
    StateSpaceDense,
)

from .ss_algorithms import (
    reduce_modal,
)


from .zpk_algorithms import (
    ss2zpk,
    zpk2cDSS,
    zpk2rDSS,
    DSS_c2r,
    zpkdict_cascade,
    poly2ss,
)


from .delay_algorithms import (
    pade_delay,
    bessel_delay,
)


from .xfer_algorithms import (
    ss2xfer,
)


from .eig_algorithms import (
    eigspaces_right,
    eigspaces_right_real,
)
