#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from .plots import IIRPlots

# REPLACE ME (if wanted)
plotter = IIRPlots()


def plot_fit(fitter, **kwargs):
    return plotter.plot_fit(fitter, **kwargs)


def plot_rel_comparison(fitter, fitter_ref, **kwargs):
    return plotter.plot_rel_comparison(fitter, fitter_ref, **kwargs)


def plot_ZP(fitter, **kwargs):
    return plotter.plot_ZP(fitter, **kwargs)


def plot_ZP(fitter, **kwargs):
    return plotter.plot_ZP(fitter, **kwargs)


def plot_ZP_S(fitter, **kwargs):
    return plotter.plot_ZP_S(fitter, **kwargs)


def plot_fitter_flag(fitter, **kwargs):
    return plotter.plot_fitter_flag(fitter, **kwargs)


def plot_fitter_flag_compare(fitter, fitter_ref, **kwargs):
    return plotter.plot_fitter_flag_compare(fitter, fitter_ref, **kwargs)


def plot_fitter_flag_residuals(fitter=None, **kwargs):
    return plotter.plot_fitter_flag_residuals(fitter, **kwargs)


def plot_ZP_grab(fitter, duals, **kwargs):
    return plotter.plot_ZP_grab(fitter, duals, **kwargs)
