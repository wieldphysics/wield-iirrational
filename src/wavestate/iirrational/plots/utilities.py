# -*- coding: utf-8 -*-
"""
"""

import numpy as np
from wavestate.bunch import Bunch


def BWz(r, F_nyquist_Hz):
    return ((1 - r) * F_nyquist_Hz)

def BWz_neg(r, F_nyquist_Hz):
    return ((r - 1) * F_nyquist_Hz)

def get_aspect(ax):
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    return (figW * w) / (figH * h)

def ZPKrep2Sf(ZPKrep):
    if ZPKrep.F_nyquist_Hz is None:
        return ZPKrep
    else:
        #TODO, make this use some other standard function
        def remap(rB):
            rB = rB.copy()
            select_r = rB.r > 0
            rB.r = (rB.r[select_r] - 1) * ZPKrep.F_nyquist_Hz
            F_Hz = np.angle(rB.c) / np.pi * ZPKrep.F_nyquist_Hz
            BW = (abs(rB.c) - 1) * ZPKrep.F_nyquist_Hz
            rB.c = BW + 1j * F_Hz
            rB.clear()
            return rB

    ret = Bunch(
        zeros         = remap(ZPKrep.zeros),
        poles         = remap(ZPKrep.poles),
        zeros_overlay = remap(ZPKrep.zeros_overlay),
        poles_overlay = remap(ZPKrep.poles_overlay),
        F_nyquist_Hz = None,
    )
    return ret

