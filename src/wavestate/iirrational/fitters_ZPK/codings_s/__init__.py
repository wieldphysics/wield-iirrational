#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import sys

import declarative
#import numpy as np
#import scipy.linalg

from ..MRF import MultiReprFilterS
from ..codings_cmn import (
    CodingGain,
    CodingDelay,
    CodingDelayNL,
)

from .cplx_real_imag import CodingRI
from .cplx_frequency_BW import CodingFBW
from .cplx_nl_freq_nl_BW import CodingnlFnlBW
#from .cplx_nl_freq_amplitude import CodingnlFA
from .real_BW import CodingRealBW
from .real_nl_BW import CodingRealnlBW
from .cplx_sos import CodingSOS
from .cplx_sos_NL import CodingSOSnl
#from .multi_order_section import CodingMOS

#get the current module object to provide a reference
_module = sys.modules[__name__]


def is_unstable_S(root):
    return root.real > 0


_common = dict(
    module = _module,
    representation = 'Sf',
    mrf_default = MultiReprFilterS,
    is_unstable = is_unstable_S,
)

coding_maps = wavestate.bunch.Bunch(
    RI = wavestate.bunch.Bunch(
        gain          = CodingGain,
        delay         = CodingDelay,
        num_r         = CodingRealBW,
        num_c         = CodingRI,
        den_r         = CodingRealBW,
        den_c         = CodingRI,
        num_r_u       = None,
        num_c_u       = None,
        den_r_u       = None,
        den_c_u       = None,
        num_collect_r = 1,
        den_collect_r = 1,
        num_collect_c = 1,
        den_collect_c = 1,
        **_common
    ),
    FBW = wavestate.bunch.Bunch(
        gain          = CodingGain,
        delay         = CodingDelay,
        num_r         = CodingRealBW,
        num_c         = CodingFBW,
        den_r         = CodingRealBW,
        den_c         = CodingFBW,
        num_r_u       = None,
        num_c_u       = None,
        den_r_u       = None,
        den_c_u       = None,
        num_collect_r = 1,
        den_collect_r = 1,
        num_collect_c = 1,
        den_collect_c = 1,
        **_common
    ),
    nlFBW = wavestate.bunch.Bunch(
        gain          = CodingGain,
        delay         = CodingDelayNL,
        num_r         = CodingRealnlBW,
        num_c         = CodingnlFnlBW,
        den_r         = CodingRealnlBW,
        den_c         = CodingnlFnlBW,
        num_r_u       = None,
        num_c_u       = None,
        den_r_u       = None,
        den_c_u       = None,
        num_collect_r = 1,
        den_collect_r = 1,
        num_collect_c = 1,
        den_collect_c = 1,
        **_common
    ),
    nlFBW_safe = wavestate.bunch.Bunch(
        gain          = CodingGain,
        delay         = CodingDelayNL,
        num_r         = CodingRealnlBW,
        num_c         = CodingnlFnlBW,
        den_r         = CodingRealnlBW,
        den_c         = CodingnlFnlBW,
        num_r_u       = CodingRealBW,
        num_c_u       = CodingFBW,
        den_r_u       = CodingRealBW,
        den_c_u       = CodingFBW,
        num_collect_r = 1,
        den_collect_r = 1,
        num_collect_c = 1,
        den_collect_c = 1,
        **_common
    ),
    #TODO, not ready, gotta be able to deal with roots at 0
    SOS = wavestate.bunch.Bunch(
        gain          = CodingGain,
        delay         = CodingDelayNL,
        num_r         = CodingSOS,
        num_c         = CodingSOS,
        den_r         = CodingSOS,
        den_c         = CodingSOS,
        num_r_u       = None,
        num_c_u       = None,
        den_r_u       = None,
        den_c_u       = None,
        num_collect_r = 2,
        den_collect_r = 2,
        num_collect_c = 1,
        den_collect_c = 1,
        **_common
    ),
    SOSnl = wavestate.bunch.Bunch(
        gain          = CodingGain,
        delay         = CodingDelayNL,
        num_r         = CodingSOSnl,
        num_c         = CodingSOSnl,
        den_r         = CodingSOSnl,
        den_c         = CodingSOSnl,
        num_r_u       = CodingSOSnl,
        num_c_u       = CodingSOSnl,
        den_r_u       = CodingSOSnl,
        den_c_u       = CodingSOSnl,
        num_collect_r = 2,
        den_collect_r = 2,
        num_collect_c = 1,
        den_collect_c = 1,
        **_common
    ),
    SOSsafe = wavestate.bunch.Bunch(
        gain          = CodingGain,
        delay         = CodingDelayNL,
        num_r         = CodingSOSnl,
        num_c         = CodingSOSnl,
        den_r         = CodingSOSnl,
        den_c         = CodingSOSnl,
        num_r_u       = CodingSOS,
        num_c_u       = CodingSOS,
        den_r_u       = CodingSOS,
        den_c_u       = CodingSOS,
        num_collect_r = 2,
        den_collect_r = 2,
        num_collect_c = 1,
        den_collect_c = 1,
        **_common
    ),
)
