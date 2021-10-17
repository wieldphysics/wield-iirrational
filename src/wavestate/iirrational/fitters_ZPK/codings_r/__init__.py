# -*- coding: utf-8 -*-
"""
"""


import declarative
#import numpy as np
#import scipy.linalg

from .base import (
    CodingType,
    Ipi,
    I2pi
)

from ..MRF import MultiReprFilterS
from .gain_delay import CodingGainDelay
from .cplx_real_imag import CodingRI
from .cplx_frequency_BW import CodingFBW
from .cplx_nl_freq_nl_BW import CodingnlFnlBW
#from .cplx_nl_freq_amplitude import CodingnlFA
from .real_BW import CodingRealBW
from .real_nl_BW import CodingRealnlBW
from .cplx_sos import CodingSOS, CodingSOSMirror
#from .multi_order_section import CodingMOS

import sys
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
        gain_delay    = CodingGainDelay,
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
        gain_delay    = CodingGainDelay,
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
        gain_delay    = CodingGainDelay,
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
        gain_delay    = CodingGainDelay,
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
        gain_delay    = CodingGainDelay,
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
    #SOS_safe = wavestate.bunch.Bunch(
    #    gain_delay    = CodingGainDelay,
    #    num_r         = CodingSOSMirror,
    #    num_c         = CodingSOSMirror,
    #    den_r         = CodingSOSMirror,
    #    den_c         = CodingSOSMirror,
    #    num_r_u       = CodingSOS,
    #    num_c_u       = CodingSOS,
    #    den_r_u       = CodingSOS,
    #    den_c_u       = CodingSOS,
    #    num_collect_r = 2,
    #    den_collect_r = 2,
    #    num_collect_c = 1,
    #    den_collect_c = 1,
    #    **_common
    #),
)
