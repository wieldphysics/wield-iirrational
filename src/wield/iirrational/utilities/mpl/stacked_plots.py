#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
from wield import bunch

from .autoniceplot import (
    asavefig,
    patch_axes,
)


def hide_xlabels(ax):
    for tl in ax.get_xticklabels():
        tl.set_visible(False)


def attach_finalizer(ax):
    ax.finalizers = []
    ax.finalized = False

    def finalize():
        ax.finalized = True
        for f in ax.finalizers:
            f()

    ax.finalize = finalize
    return


def generate_ax(ax_group=None):
    ax = bunch.DeepBunch()
    attach_finalizer(ax)
    if ax_group is not None:
        ax_group.finalizers.append(ax.finalize)
    return ax


def generate_stacked_plot_ax(
    name_use_list,
    gs_base=gridspec.GridSpec(1, 1)[0],
    heights_phys_in_default=1.5,
    heights_phys_in={},
    height_ratios={},
    width_ratios=[1],
    xscales="linear",
    xlim=None,
    wspacing=0.04,
    width_phys_in=None,
    fig=None,
    ax_group=None,
    hspace=0.2,
    **kwargs
):
    view_names = []
    height_ratio_list = []
    height_phys_in = 0
    autocalls = dict()
    for ntup in name_use_list:
        if isinstance(ntup, (tuple, list)):
            name, b_use = ntup
        else:
            name = ntup
            b_use = True
        if callable(name):
            autocalls[name] = b_use
            name = name.__name__
        else:
            if not b_use:
                continue
        height_ratio = height_ratios.get(name, 1)
        height_ratio_list.append(height_ratio)
        height_phys_in += height_ratio * heights_phys_in.get(
            name, heights_phys_in_default
        )
        view_names.append(name)

    gs_DC = gridspec.GridSpecFromSubplotSpec(
        len(view_names),
        len(width_ratios),
        subplot_spec=gs_base,
        height_ratios=height_ratio_list,
        width_ratios=width_ratios,
        wspace=wspacing,
        hspace=hspace,
    )

    if not isinstance(xscales, (list, tuple)):
        xscales = [xscales] * len(width_ratios)

    if width_phys_in is None:
        width_phys_in = mpl.rcParams["figure.figsize"][0]

    if fig is None:
        fig = plt.figure()
        fig.set_size_inches(width_phys_in, height_phys_in)
        # fig.set_dpi(160)

    def hide_finalizer(axB, ax_local):
        def finalize():
            hide_xlabels(ax_local)

        return finalize

    def scale_finalizer(ax_local, scale):
        def finalize():
            ax_local.set_xscale(scale)

        return finalize

    axB = generate_ax(ax_group)
    axB.fig = fig

    def save(rootname, **kwargs):
        axB.finalize()
        axB << asavefig(rootname, **kwargs)

    axB.save = save

    for col_idx in range(len(width_ratios)):
        ax_top = None
        ax_list = []
        for idx, name in enumerate(view_names):
            ax_local = fig.add_subplot(gs_DC[idx, col_idx], sharex=ax_top)
            # patch_axes(ax_local)
            ax_local.grid(b=True)
            ax_local.grid(b=True, which="minor", color=(0.9, 0.9, 0.9), lw=0.5)
            axB["ax{0}_{1}".format(idx, col_idx)] = ax_local
            if col_idx == 0:
                axB["ax{0}".format(idx)] = ax_local
            axB[name + "_{0}".format(col_idx)] = ax_local
            if col_idx == 0:
                axB[name] = ax_local
            if ax_top is None:
                ax_top = ax_local
            ax_list.append(ax_local)
            axB.finalizers.append(scale_finalizer(ax_local, xscales[col_idx]))
            if idx < len(view_names) - 1:
                axB.finalizers.append(hide_finalizer(axB, ax_local))
        axB["ax_bottom_{0}".format(col_idx)] = ax_local
        axB["ax_top_{0}".format(col_idx)] = ax_top
        axB["ax_list_{0}".format(col_idx)] = ax_list
        if col_idx == 0:
            axB["ax_bottom"] = ax_local
            axB["ax_top"] = ax_top
            axB["ax_list"] = ax_list
    for call, cparams in autocalls.items():
        call(axB[call.__name__], **cparams)
    return axB
