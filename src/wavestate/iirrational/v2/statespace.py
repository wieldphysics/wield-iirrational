# -*- coding: utf-8 -*-
"""
"""

import numpy as np

from ..TFmath.statespace import (
    ss2zpk,
    ss2xfer
)
from .. import fitters_ZPK
from . import data2filter


def ss2filter(
    A, B, C, D,
    F_Hz,
    SNR = 1e3,
    idx_in = None,
    idx_out = None,
    mode = 'reduce',
    prune_Qrank = .01,
    prune_Qrank_unstable = .5,
    **kwargs
):
    z, p, k = ss2zpk(
        A, B, C, D,
        idx_in = idx_in,
        idx_out = idx_out,
        fmt = 'IIRrational',
        Q_rank_cutoff_unstable = prune_Qrank_unstable,
    )

    fit = data2filter(
        data = ss2xfer(A, B, C, D, F_Hz, idx_in = idx_in, idx_out = idx_out),
        F_Hz = F_Hz,
        zeros = z,
        poles = p,
        gain  = k,
        SNR = SNR,
        mode = mode,
        prune_Qrank = prune_Qrank,
        coding_map = fitters_ZPK.coding_maps.RI,
        **kwargs
    )
    return fit
