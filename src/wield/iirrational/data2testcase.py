#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from wield.bunch import Bunch
from .utilities import args
import scipy.io


def data2testcase(
    argB=args.UNSPEC,
    data=args.UNSPEC,
    F_Hz=args.UNSPEC,
    SNR=args.UNSPEC,
    fname=args.UNSPEC,
    F_nyquist_Hz=args.UNSPEC,
    bestfit_ZPK_z=args.UNSPEC,
    bestfit_ZPK_s=args.UNSPEC,
    description=args.UNSPEC,
):
    """
    Save a testcase in in matlab .mat format with a standardized .dat for standard read testcase reading.
    :arg: data
    """
    data = args.argscan(locals(), argB, args.REQ, arg="data")
    F_Hz = args.argscan(locals(), argB, args.REQ, arg="F_Hz")
    SNR = args.argscan(locals(), argB, 1, arg="SNR")
    fname = args.argscan(locals(), argB, args.REQ, arg="fname")
    F_nyquist_Hz = args.argscan(locals(), argB, args.REQ, arg="F_nyquist_Hz")
    bestfit_ZPK_z = args.argscan(locals(), argB, None, arg="bestfit_ZPK_z")
    bestfit_ZPK_s = args.argscan(locals(), argB, None, arg="bestfit_ZPK_s")
    description = args.argscan(locals(), argB, args.REQ, arg="description")

    mdict = dict(
        data=data,
        F_Hz=F_Hz,
        SNR=SNR,
        description=description,
        version=1,
    )

    if bestfit_ZPK_s is not None:
        mdict["bestfit_ZPK_s"] = bestfit_ZPK_s

    if bestfit_ZPK_z is not None:
        mdict["bestfit_ZPK_z"] = bestfit_ZPK_z

    return scipy.io.savemat(fname, mdict)


def testcase2data(fname):
    import scipy.io

    data = scipy.io.loadmat(fname, squeeze_me=True)
    return Bunch(data)
