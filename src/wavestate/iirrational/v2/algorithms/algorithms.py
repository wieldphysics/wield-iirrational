# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np


from .MRFprogram import (
    resrank_program,
    ranking_reduction_trials,
)

def sign_check_flip(fitter, max_N = None):
    """
    """
    if max_N is None:
        xfer = fitter.xfer_fit
        data = fitter.data
        W = fitter.W
    else:
        xfer = fitter.xfer_fit[:max_N]
        data = fitter.data[:max_N]
        W = fitter.W[:max_N]
    rat = data / xfer
    select = np.isfinite(rat)
    rat_ang = np.exp(1j * np.angle(rat))
    ang_avg_fit = np.sum(rat_ang[select] * W[select]**2) / np.sum(W[select]**2)

    residuals = fitter.residuals_average
    if (ang_avg_fit.real < 0):
        fitter.gain = -fitter.gain

    if (fitter.residuals_average > residuals):
        fitter.gain = -fitter.gain

    if max_N is None:
        xfer = fitter.xfer_fit
        data = fitter.data
        W = fitter.W
    else:
        xfer = fitter.xfer_fit[:max_N]
        data = fitter.data[:max_N]
        W = fitter.W[:max_N]
    rat = data / xfer
    rat_ang = np.exp(1j * np.angle(rat))
    ang_avg_fit = np.sum(rat_ang[select] * W[select]**2) / np.sum(W[select]**2)
    assert(ang_avg_fit.real > 0)
    return


def optimize_anneal(aid, fitter = None):
    if fitter is None:
        fitter = aid.fitter
    try:
        #TODO, maybe try nelder_mead here
        fitter.optimize(
            residuals_type = aid.hint('residuals_type_alt', 'log'),
            aid            = aid,
            max_nfev       = 30,
        )
    except Exception as e:
        aid.log_debug(10, "Optimize Exception, ", e)

    try:
        fitter.optimize(aid = aid)
    except Exception as e:
        aid.log_debug(10, "Optimize Exception, ", e)
        return False
    return True



