# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
import scipy.signal


def SNR_estimate(F_Hz, xfer, width = 5):
    #xfer = xfer / abs(xfer)

    xferL = np.log(abs(xfer))
    d = np.cumsum(xferL)
    diff = (d[width:] - d[:-width])/width
    diff = np.concatenate([diff, d[-1] - d[-width-1:-1] / np.arange(1 + width, 1, -1), ])
    R1 = xferL - diff
    R1d = R1[1:] - R1[:-1]

    xferP = xfer/abs(xfer)
    d = np.cumsum(xferP)
    diff = (d[width:] - d[:-width])/width
    diff = np.concatenate([diff, d[-1] - d[-width-1:-1] / np.arange(1 + width, 1, -1), ])
    R2 = xferP / diff
    R2d = R2[1:] - R2[:-1]

    W_c = np.cumsum(np.maximum(abs(R2d)**2, abs(R1d)**2))
    W_diff = (W_c[width:] - W_c[:-width])/width
    W_est = np.maximum(1/W_diff - 0.5, 0)**0.5
    W_est = np.concatenate([W_est, [W_est[-1]] * (width + 1), ])

    return W_est


