#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from matplotlib import gridspec
import matplotlib.pyplot as plt

from .. import fitters_ZPK

from ..utilities.mpl import (
    asavefig,
    generate_ax,
)


def plot_fitter_flag(
    self,
    fitter,
    title=None,
    plot_zp=True,
    plot_to_nyquist=True,
    with_error=True,
    with_data_error=True,
    N_points=300,
    gs_base=gridspec.GridSpec(1, 1)[0],
    fig=None,
    ax_bunch=None,
    xscale=None,
    plot_data_scale=10,
    plot_data_scale_low=2,
):
    if not isinstance(fitter, fitters_ZPK.MultiReprFilterBase):
        fitter = fitters_ZPK.ZPKrep2MRF(
            fitter.ZPKrep,
            coding_map=fitters_ZPK.coding_maps.RI,
        )
    axB = generate_ax(ax_bunch)

    if fitter.F_nyquist_Hz is None or fitter.F_nyquist_Hz > 1.2 * fitter.F_max_Hz:
        if fig is None:
            fig = plt.figure()
            fig.set_size_inches(12, 9)
            axB.fig = fig
        gs_DC = gridspec.GridSpecFromSubplotSpec(
            2,
            2,
            subplot_spec=gs_base,
            wspace=0.2,
            hspace=0.2,
            # height_ratios = height_ratio_list,
            width_ratios=[1, 0.3],
        )

        gs_data = gs_DC[0, 0]
        gs_nyquist = gs_DC[1, 0]
        gs_zp = gs_DC[:, 1]
        zp_horz = False
    else:
        if fig is None:
            fig = plt.figure()
            fig.set_size_inches(10, 9)
            axB.fig = fig
        gs_DC = gridspec.GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=gs_base,
            wspace=0.1,
            hspace=0.3,
            height_ratios=[1, 0.7],
        )

        gs_data = gs_DC[0, 0]
        gs_nyquist = False
        gs_zp = gs_DC[1, 0]
        zp_horz = True

    if title is not None:
        axB.fig.suptitle(title)

    if gs_data:
        axB_top = self.plot_fit(
            fitter=fitter,
            ax_bunch=None,
            gs_base=gs_data,
            with_error=with_error,
            with_data_error=with_data_error,
            plot_to_nyquist=False,
            plot_now=False,
            xscale=xscale,
            fig=fig,
        )
        axB.top = axB_top
        axB_top.plot_data(fitter)
        axB_top.plot_fit(fitter)
        axB_top.plot_zp(fitter)
        axB.finalizers.append(axB.top.finalize)

    if gs_nyquist:
        axB_bottom = self.plot_fit(
            fitter=fitter,
            ax_bunch=None,
            gs_base=gs_nyquist,
            with_error=False,
            with_data_error=False,
            plot_past_data=True,
            plot_to_nyquist=True,
            plot_now=False,
            fig=fig,
            xscale="log",
            plot_data_scale=plot_data_scale,
            plot_data_scale_low=plot_data_scale_low,
        )
        axB.bottom = axB_bottom
        axB_bottom.plot_data(fitter)
        axB_bottom.plot_fit(fitter)
        axB_bottom.plot_zp(fitter)
        left, right = axB.bottom.ax_bottom.get_xlim()
        axB.finalizers.append(axB.bottom.finalize)
        if fitter.F_nyquist_Hz is not None:
            axB.bottom.ax_bottom.set_xlim(left, fitter.F_nyquist_Hz)

    if gs_zp:
        if fitter.F_nyquist_Hz is None or fitter.F_nyquist_Hz / max(fitter.F_Hz) > 5:
            axB_right = self.plot_ZP_S(
                fitter=fitter,
                ax_bunch=None,
                gs_base=gs_zp,
                horizontal=zp_horz,
                fig=fig,
            )
        else:
            axB_right = self.plot_ZP(
                fitter=fitter,
                ax_bunch=None,
                gs_base=gs_zp,
                horizontal=zp_horz,
                fig=fig,
            )
        axB.right = axB_right

    def plot_fit(*args, **kwargs):
        db = axB.top.plot_fit(*args, **kwargs)
        axB.bottom.plot_fit(*args, db_ref=db, **kwargs)

    axB.plot_fit = plot_fit

    def save(rootname, **kwargs):
        axB.finalize()
        axB << asavefig(rootname, **kwargs)

    axB.save = save
    axB.finalize()
    return axB


def plot_fitter_flag_compare(
    self,
    fitter,
    fitter_ref,
    N_points=300,
    gs_base=gridspec.GridSpec(1, 1)[0],
    fig=None,
    ax_bunch=None,
):
    if not isinstance(fitter, fitters_ZPK.MultiReprFilterBase):
        fitter = fitters_ZPK.RDF2MRF(fitter)
    if not isinstance(fitter_ref, fitters_ZPK.MultiReprFilterBase):
        fitter_ref = fitters_ZPK.RDF2MRF(fitter_ref)
    axB = generate_ax(ax_bunch)

    if fig is None:
        fig = plt.figure()
        fig.set_size_inches(14, 9)
        axB.fig = fig
    gs_DC = gridspec.GridSpecFromSubplotSpec(
        3,
        3,
        subplot_spec=gs_base,
        wspace=0.3,
        hspace=0.2,
        height_ratios=[1, 1, 0.5],
        width_ratios=[1, 0.3, 0.3],
    )

    gs_data = gs_DC[0, 0]
    gs_nyquist = gs_DC[1, 0]
    gs_zp = gs_DC[:2, 1]
    gs_zp_ref = gs_DC[:2, 2]
    zp_horz = False
    gs_rel = gs_DC[2, :]

    if gs_data:
        axB_top = self.plot_fit(
            fitter=fitter,
            gs_base=gs_data,
            with_error=True,
            with_data_error=True,
            plot_to_nyquist=False,
            fig=fig,
            ax_bunch=None,
        )
        axB.top = axB_top
        axB_top = self.plot_fit(
            fitter=fitter_ref,
            gs_base=gs_data,
            with_error=False,
            with_data_error=False,
            plot_data=False,
            plot_to_nyquist=False,
            fig=fig,
            ax_bunch=axB_top,
            plot_zp=False,
            label_fit="original ({nP}P, {nZ}Z, Rsq={Rsq:.2e})",
        )
        axB_top.mag.legend(ncol=1, fontsize=8, loc="best").set_visible(False)

    if gs_nyquist:
        axB_bottom = self.plot_fit(
            fitter=fitter,
            ax_bunch=None,
            gs_base=gs_nyquist,
            with_error=False,
            with_data_error=False,
            plot_to_nyquist=True,
            fig=fig,
            xscale="log",
        )
        axB_bottom = self.plot_fit(
            fitter=fitter_ref,
            gs_base=gs_nyquist,
            with_error=False,
            with_data_error=False,
            plot_data=False,
            plot_to_nyquist=True,
            fig=fig,
            xscale="log",
            ax_bunch=axB_bottom,
            plot_zp=False,
            label_fit="original ({nP}P, {nZ}Z, Rsq={Rsq:.2e})",
        )
        axB.bottom = axB_bottom

        def finalize():
            axB_bottom.mag.legend(ncol=1, fontsize=8, loc="best")
            left, right = axB.bottom.ax_bottom.get_xlim()
            axB.bottom.ax_bottom.set_xlim(left, fitter.F_nyquist_Hz)

        axB_bottom.finalizers.append(finalize)

    if gs_zp:
        if fitter.F_nyquist_Hz / max(fitter.F_Hz) > 5:
            axB_right = self.plot_ZP_S(
                fitter=fitter,
                ax_bunch=None,
                gs_base=gs_zp,
                horizontal=zp_horz,
                fig=fig,
            )
        else:
            axB_right = self.plot_ZP(
                fitter=fitter,
                ax_bunch=None,
                gs_base=gs_zp,
                horizontal=zp_horz,
                fig=fig,
            )
        axB.right = axB_right
        axB.right.ax_poles.set_title("fit roots (s)")

    if gs_zp_ref:
        if fitter.F_nyquist_Hz / max(fitter.F_Hz) > 5:
            axB_right = self.plot_ZP_S(
                fitter=fitter_ref,
                ax_bunch=None,
                gs_base=gs_zp_ref,
                horizontal=zp_horz,
                fig=fig,
            )
        else:
            axB_right = self.plot_ZP(
                fitter=fitter_ref,
                ax_bunch=None,
                gs_base=gs_zp_ref,
                horizontal=zp_horz,
                fig=fig,
            )
        axB.rright = axB_right
        axB.rright.ax_poles.set_title("reference roots (s)")

    if gs_rel:
        self.plot_rel_comparison(
            fitter=fitter,
            fitter_ref=fitter_ref,
            gs_base=gs_rel,
            fig=fig,
        )

    def save(rootname, **kwargs):
        axB << asavefig(rootname, **kwargs)

    axB.save = save

    axB.finalize()
    return axB


def plot_fitter_flag_residuals(
    self,
    fitter,
    title=None,
    plot_zp=True,
    plot_to_nyquist=True,
    with_error=True,
    with_data_error=True,
    N_points=300,
    gs_base=gridspec.GridSpec(1, 1)[0],
    fig=None,
    ax_bunch=None,
    xscale=None,
):
    if not isinstance(fitter, fitters_ZPK.MultiReprFilterBase):
        fitter = fitters_ZPK.ZPKrep2MRF(
            fitter.ZPKrep,
            coding_map=fitters_ZPK.coding_maps.RI,
        )
    axB = generate_ax(ax_bunch)

    if fig is None:
        fig = plt.figure()
        fig.set_size_inches(12, 9)
        axB.fig = fig

    gs_DC = gridspec.GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=gs_base,
        wspace=0.2,
        hspace=0.2,
    )

    gs_data = gs_DC[0, 0]
    gs_res = gs_DC[1, 0]

    if title is not None:
        axB.fig.suptitle(title)

    if gs_data:
        axB_top = self.plot_fit(
            fitter=fitter,
            ax_bunch=None,
            gs_base=gs_data,
            with_error=with_error,
            with_data_error=with_data_error,
            plot_past_data=True,
            plot_to_nyquist=True,
            plot_now=False,
            fig=fig,
            xscale=xscale,
        )
        axB_top.plot_data(fitter)
        db_fit = axB_top.plot_fit(fitter)
        axB_top.plot_zp(fitter)
        axB.top = axB_top
        axB.finalizers.append(axB.top.finalize)

    if gs_res:
        #axB_bottom = self.plot_residuals(
        #    fitter=fitter,
        #    ax_bunch=None,
        #    gs_base=gs_res,
        #    plot_now=False,
        #    fig=fig,
        #    xscale=xscale,
        #)
        axB_bottom = self.plot_fit(
            fitter=fitter,
            ax_bunch=None,
            gs_base=gs_res,
            with_error=with_error,
            with_data_error=with_data_error,
            plot_past_data=False,
            plot_to_nyquist=True,
            plot_now=False,
            residuals=True,
            fig=fig,
            xscale=xscale,
        )
        axB.bottom = axB_bottom
        axB_bottom.mag.axhline(1.1, ls='--', color='black', lw=1)
        axB_bottom.mag.axhline(1.01, ls='--', color='black', lw=1)
        axB_bottom.mag.axhline(1/1.1, ls='--', color='black', lw=1)
        axB_bottom.mag.axhline(1/1.01, ls='--', color='black', lw=1)
        axB_bottom.plot_data(fitter)
        axB_bottom.plot_fit(fitter, db_ref=db_fit)
        left, right = axB.bottom.ax_bottom.get_xlim()
        if fitter.F_nyquist_Hz is not None:
            axB.bottom.ax_bottom.set_xlim(left, fitter.F_nyquist_Hz)
        axB.finalizers.append(axB.bottom.finalize)

    def plot_fit(*args, **kwargs):
        db = axB.top.plot_fit(*args, **kwargs)
        axB.bottom.plot_fit(*args, db_ref=db, **kwargs)

    axB.plot_fit = plot_fit

    def save(rootname, **kwargs):
        axB.finalize()
        axB << asavefig(rootname, **kwargs)

    axB.save = save
    axB.finalize()
    return axB
