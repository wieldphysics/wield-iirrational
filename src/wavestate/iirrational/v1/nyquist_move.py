#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import numpy as np
import declarative

from ..utilities import ensure_aid
from .. import TFmath
from .. import representations


def nyquist_move(
    fitter,
    nyquist_new,
    clear_neg_real = False,
    split_neg_real = False,
    aid            = None,
):
    """
    Modifies a fitter IN PLACE to use a different nyquist frequency
    """
    #TODO, update this to properly use RootBunches and ZPKreps
    aid = ensure_aid(aid)
    ZPKrep = fitter.ZPKrep

    F_nyquist_Hz = ZPKrep.F_nyquist_Hz
    op = wavestate.bunch.Bunch(
        gain_flip = 1
    )

    #currently a restriction in the algorithm, it can't handle aliasing of poles (need bilinear for that)
    assert(nyquist_new > F_nyquist_Hz)
    ratio = nyquist_new / F_nyquist_Hz

    def root_list_map(root_list, rtype):
        new_list = []
        for root in root_list:
            new_roots = None
            if abs(root.imag) < 1e-15:
                #aid.log(root)
                if root.real > 0:
                    if root.real < 1:
                        new_roots = [1 - (1 - root.real) / ratio]
                    else:
                        #actually, just clear it
                        aid.log("Cleared ", rtype, root)
                        #new_roots = []
                        #flip it then insert
                        op.gain_flip *= -1
                        root = 1/root
                        #new_roots = [1 - (1 - root.real) / ratio]
                        new_roots = []
                else:
                    if clear_neg_real:
                        aid.log("Cleared ", rtype, root)
                        new_roots = []
                    elif split_neg_real:
                        angle = np.pi / ratio
                        amp = -root.real
                        if amp < 1:
                            amp = 1 - (1 - amp) / ratio
                        else:
                            amp = (amp - 1) / ratio + 1
                        new_roots = [
                            amp * np.exp(1j * angle),
                            amp * np.exp(-1j * angle)
                        ]
                    else:
                        new_roots = [-1 - (-1 - root.real) / ratio]
            else:
                amp = abs(root)
                angle = np.angle(root)
                if angle > 0:
                    angle = angle / ratio
                else:
                    angle = angle / ratio
                if amp < 1:
                    amp = 1 - (1 - amp) / ratio
                else:
                    amp = (amp - 1) / ratio + 1
                new_roots = [amp * np.exp(1j * angle)]
            new_list.extend(new_roots)
        return new_list

    xfer_pre = fitter.xfer_fit

    nyquist_prev = ZPKrep.F_nyquist_Hz
    fitter.F_nyquist_Hz = nyquist_new
    poles = root_list_map(fitter.poles.fullplane, rtype = 'pole')
    zeros = root_list_map(fitter.zeros.fullplane, rtype = 'zero')
    fitter.poles = poles
    fitter.zeros = zeros

    #don't use the root mapper above since it does special things and can drop
    #overlays
    Zov, Pov, Kex = TFmath.SorZtoSorZ(
        (
            fitter.zeros_overlay.fullplane,
            fitter.poles_overlay.fullplane,
            1
        ),
        F_nyquist_Hz_in  = nyquist_prev,
        F_nyquist_Hz_out = nyquist_new,
    )
    fitter.poles_overlay = Pov
    fitter.zeros_overlay = Zov

    xfer_post = fitter.xfer_fit

    gain_scale = np.median(abs(xfer_post / xfer_pre))

    gain = fitter.gain / gain_scale * op.gain_flip
    fitter.gain = gain
    #aid.log("GAIN: ", fitter.gain)
    return fitter


def nyquist_remove(
    fitter,
    clear_neg_real = False,
    split_neg_real = False,
    aid            = None,
):
    """
    Modifies a fitter IN PLACE to use a different nyquist frequency
    """
    aid = ensure_aid(aid)

    F_nyquist_Hz = fitter.F_nyquist_Hz
    op = wavestate.bunch.Bunch(
        gain_flip = 1
    )

    def root_list_map(root_list, rtype):
        new_list = []
        for root in root_list:
            new_roots = None
            if abs(root.imag) < 1e-15:
                #aid.log(root)
                if root.real > 0:
                    if root.real < 1:
                        new_roots = [F_nyquist_Hz * (root.real - 1)]
                    else:
                        aid.log("Cleared ", rtype, root)
                        op.gain_flip *= -1
                        new_roots = []
                else:
                    if clear_neg_real:
                        aid.log("Cleared ", rtype, root)
                        #TODO, not sure about this one
                        op.gain_flip *= -1
                        new_roots = []
                    elif split_neg_real:
                        BW = -(1 + root.real) * F_nyquist_Hz
                        new_roots = [
                            BW + 1j * F_nyquist_Hz,
                            BW + -1j * F_nyquist_Hz,
                        ]
                    else:
                        new_roots = [F_nyquist_Hz * -(root.real + 1)]
            else:
                amp = abs(root)
                F_Hz = np.angle(root) / np.pi * F_nyquist_Hz
                new_roots = [(amp - 1) * F_nyquist_Hz + 1j * F_Hz]
            new_list.extend(new_roots)
        return new_list

    xfer_pre = fitter.xfer_fit

    poles = root_list_map(fitter.poles.fullplane, rtype = 'pole')
    zeros = root_list_map(fitter.zeros.fullplane, rtype = 'zero')

    #don't use the root mapper above since it does special things and can drop
    #overlays
    Zov, Pov, Kex = TFmath.ZtoS(
        (
            fitter.zeros_overlay.fullplane,
            fitter.poles_overlay.fullplane,
            1
        ),
        F_nyquist_Hz = fitter.F_nyquist_Hz,
    )

    ZPKrep = representations.ZPKwData(
        parent        = fitter.ZPKrep,
        poles         = poles,
        zeros         = zeros,
        poles_overlay = Zov,
        zeros_overlay = Pov,
        gain          = Kex * fitter.ZPKrep.gain
    )

    xfer_post = ZPKrep.xfer_fit

    gain_scale = np.median(abs(xfer_post / xfer_pre))
    gain = fitter.gain / gain_scale * op.gain_flip

    ZPKrep.gain = gain
    return ZPKrep
