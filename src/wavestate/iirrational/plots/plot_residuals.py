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
from matplotlib import gridspec

from ..utilities.mpl import generate_stacked_plot_ax
from .plot_augment import plot_augment_residuals


def plot_residuals(
    self,
    fitter     = None,
    gs_base    = gridspec.GridSpec(1, 1)[0],
    fig        = None,
    ax_bunch   = None,
    xscale     = None,
    plot_now   = True,
    **kwargs
):
    use_phase = True

    if ax_bunch is None:
        axB = generate_stacked_plot_ax(
            [
                ('mag', True),
                ('phase', use_phase),
            ],
            height_ratios = dict(
                mag = 2,
            ),
            heights_phys_in_default = 1,
            xscales = 'log',
            gs_base = gs_base,
            fig     = fig,
            ax_bunch = ax_bunch,
        )
    else:
        axB = ax_bunch

    fit_list = []
    axB.plot_fit = plot_augment_residuals(
        self,
        ax_mag = axB.mag,
        ax_phase = axB.phase,
        fit_list = fit_list,
        **kwargs
    )

    axB.mag.set_ylabel('Residuals Magnitude')
    if axB.phase:
        axB.phase.set_ylabel('Residuals Phase [deg]')
    axB.ax_bottom.set_xlabel('Frequency [Hz]')

    if fitter is not None and plot_now:
        axB.plot_fit(fitter)

    def finalize():
        min_lims = []
        max_lims = []
        for fitB in fit_list:
            min_err_lim = min(
                fitB.median_res_lim,
                fitB.min_res_err_lim,
                max(fitB.min_res_err_lim / 30, fitB.min_res_err_lim * 1)
            )
            min_lims.append(min_err_lim)
            max_err_lim = max(
                fitB.median_res_lim,
                fitB.max_res_err_lim,
                min(fitB.max_res_err_lim * 30, fitB.max_res_err_lim / 1)
            )
            max_lims.append(max_err_lim)
        axB.mag.set_ylim(np.min(min_lims), max(max_lims))

        scales = []
        mins   = []
        minsNZ = []
        maxs   = []
        for fitB in fit_list:
            scales.append(fitB.scale_log)
            mins.append(fitB.minF_Hz)
            minsNZ.append(fitB.minF_HzNZ)
            maxs.append(fitB.maxF_Hz)

        if xscale is None:
            if np.count_nonzero(scales) / len(scales) > .5:
                final_xscale = 'log_zoom'
            else:
                final_xscale = 'linear'
        else:
            final_xscale = xscale

        if final_xscale in ['log', 'log_zoom']:
            xmin = np.min(minsNZ)
        else:
            xmin = np.min(mins)
        xmax = np.max(maxs)
        axB.ax_bottom.set_xscale(final_xscale)
        axB.ax_bottom.set_xlim(xmin, xmax)
        for subax in axB.ax_list_0:
            subax.set_xscale(final_xscale)
            subax.set_xlim(xmin, xmax)
        axB.mag.legend(
            ncol = 2,
            fontsize = 8,
            loc = 'best'
        )

    axB.finalizers.append(finalize)
    return axB
