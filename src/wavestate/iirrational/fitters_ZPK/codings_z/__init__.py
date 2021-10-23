#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


# import numpy as np
# import scipy.linalg
from wavestate import declarative

from ..codings_cmn import (
    CodingTypeZ,
    Ipi,
    I2pi,
    BranchCutAmbiguity,
    CodingGain,
    CodingDelay,
)

from ..MRF import MultiReprFilterZ
from .cplx_real_imag import CodingRI

from .cplx_frequency_BW import CodingFBW
from .cplx_freq_nl_BW import CodingFnlBW
from .cplx_nl_freq_nl_BW import CodingnlFnlBW

from .cplx_frequency_amplitude import CodingFA, CodingFAUX
from .cplx_frequency_nl_amp import CodingFnlA
from .cplx_nl_freq_nl_amp import CodingnlFnlA

from .real_amplitude import CodingRealAmp, CodingRealAmpUX
from .real_BW import CodingRealBW

from .real_nl_amplitude import CodingRealnlA
from .real_nl_BW import CodingRealnlBW

from .cplx_sos import CodingSOS, CodingSOSMirror
from .multi_order_section import CodingMOS

import sys

# get the current module object to provide a reference
_module = sys.modules[__name__]


def is_unstable_Z(root):
    return abs(root) >= 1


_common = dict(
    module=_module,
    representation="Z",
    mrf_default=MultiReprFilterZ,
    is_unstable=is_unstable_Z,
)


coding_maps = wavestate.bunch.Bunch(
    nlFBW=wavestate.bunch.Bunch(
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingRealnlBW,
        num_c=CodingnlFnlBW,
        den_r=CodingRealnlBW,
        den_c=CodingnlFnlBW,
        num_r_u=None,
        num_c_u=None,
        den_r_u=None,
        den_c_u=None,
        num_collect_r=1,
        den_collect_r=1,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
    nlFBW_safe=wavestate.bunch.Bunch(
        # gain    = CodingGainNLDelay,
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingRealnlBW,
        num_c=CodingnlFnlBW,
        den_r=CodingRealnlBW,
        den_c=CodingnlFnlBW,
        num_r_u=CodingRealBW,
        num_c_u=CodingFBW,
        den_r_u=CodingRealBW,
        den_c_u=CodingFBW,
        num_collect_r=1,
        den_collect_r=1,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
    nlFA=wavestate.bunch.Bunch(
        # gain    = CodingGainNLDelay,
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingRealnlA,
        num_c=CodingnlFnlA,
        den_r=CodingRealnlA,
        den_c=CodingnlFnlA,
        num_r_u=None,
        num_c_u=None,
        den_r_u=None,
        den_c_u=None,
        num_collect_r=1,
        den_collect_r=1,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
    nlFA_safe=wavestate.bunch.Bunch(
        # gain    = CodingGainNLDelay,
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingRealnlA,
        num_c=CodingnlFnlA,
        den_r=CodingRealnlA,
        den_c=CodingnlFnlA,
        num_r_u=CodingRealAmp,
        num_c_u=CodingFA,
        den_r_u=CodingRealAmp,
        den_c_u=CodingFA,
        num_collect_r=1,
        den_collect_r=1,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
    FA=wavestate.bunch.Bunch(
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingRealAmp,
        num_c=CodingFA,
        den_r=CodingRealAmp,
        den_c=CodingFA,
        num_r_u=None,
        num_c_u=None,
        den_r_u=None,
        den_c_u=None,
        num_collect_r=1,
        den_collect_r=1,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
    FBW=wavestate.bunch.Bunch(
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingRealBW,
        num_c=CodingFBW,
        den_r=CodingRealBW,
        den_c=CodingFBW,
        num_r_u=None,
        num_c_u=None,
        den_r_u=None,
        den_c_u=None,
        num_collect_r=1,
        den_collect_r=1,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
    RI=wavestate.bunch.Bunch(
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingRealAmp,
        num_c=CodingRI,
        den_r=CodingRealAmp,
        den_c=CodingRI,
        num_r_u=None,
        num_c_u=None,
        den_r_u=None,
        den_c_u=None,
        num_collect_r=1,
        den_collect_r=1,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
    SOS=wavestate.bunch.Bunch(
        gain=CodingGain,
        delay=CodingDelay,
        num_r=CodingSOS,
        num_c=CodingSOS,
        den_r=CodingSOS,
        den_c=CodingSOS,
        num_r_u=None,
        num_c_u=None,
        den_r_u=None,
        den_c_u=None,
        num_collect_r=2,
        den_collect_r=2,
        num_collect_c=1,
        den_collect_c=1,
        **_common
    ),
)
