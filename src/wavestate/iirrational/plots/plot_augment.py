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
from wavestate.bunch import Bunch

from ..utilities import args

from .utilities import ZPKrep2Sf
from ..utilities.np import domain_sort
from .. import fitters_ZPK


def angle_cutter(
    F_Hz,
    data,
    deg=False,
    semiroll=False,
):
    F_Hz, data = domain_sort(F_Hz, data)
    ang = np.angle(data)

    above = ang > np.pi / 2
    below = ang < -np.pi / 2
    jumps = (above[:-1] & below[1:]) | (below[:-1] & above[1:])

    jump_locs = np.nonzero(jumps)[0]

    if semiroll:
        for idx1, idx2 in zip(jump_locs[:-1] + 1, jump_locs[1:]):
            if np.all(above[idx1:idx2]):
                ang[idx1:idx2] -= np.pi
            elif np.all(below[idx1:idx2]):
                ang[idx1:idx2] += np.pi
            else:
                # it is crossing, do nothing
                pass

    ins_locs = jump_locs + 1
    ins_F_Hz = (F_Hz[jump_locs] + F_Hz[jump_locs]) / 2
    F_Hz = np.insert(F_Hz, ins_locs, ins_F_Hz)
    if deg:
        ang = ang * 180 / np.pi
    ang = np.insert(ang, ins_locs, np.nan)
    return (F_Hz, ang)


def plot_augment_data(
    self,
    ax_mag=None,
    ax_phase=None,
    data_list=[],
    with_error=True,
    label=args.UNSPEC,
    markersize=3,
    marker=".",
    linestyle="",
    **kwargs
):
    label_default = args.argscan(label, self.label_data, "data")

    def plot(fitter, label=None, db_ref=None, color=None, **kwargs):
        if label is None:
            label = label_default
        db = Bunch()
        db.rep = fitter

        if color is None and db_ref is not None:
            color = db_ref.mline.get_color()

        if marker is not None:
            kwargs["marker"] = marker
        if markersize is not None:
            kwargs["markersize"] = markersize
        if linestyle is not None:
            kwargs["linestyle"] = linestyle

        f = fitter.F_Hz
        h = fitter.data
        if with_error:
            W = 0.001 + fitter.W
            select = f > 0
            db.min_data_err_lim = np.min(abs(h)) * max(1 / 3, np.min(1 / (1 + 1 / W)))
            db.max_data_err_lim = np.max(abs(h)) * min(3, np.max((1 + 1 / W)))
            if ax_mag:
                # ax_mag.fill_between(
                #    fitter.F_Hz[select],
                #    np.maximum((abs(h) / (1 + 1 / W))[select], np.min(abs(h)) / 100),
                #    np.minimum((abs(h) * (1 + 1 / W))[select], np.max(abs(h)) * 100),
                #    color = color,
                #    alpha = .2,
                #    **kwargs
                # )
                abs_h = abs(h)
                ax_mag.errorbar(
                    fitter.F_Hz[select],
                    abs_h[select],
                    yerr=(
                        abs_h[select]
                        - np.maximum(
                            (abs_h / (1 + 1 / W))[select], np.min(abs_h[select]) / 100
                        ),
                        np.minimum(
                            (abs_h * (1 + 1 / W))[select], np.max(abs_h[select]) * 100
                        )
                        - abs_h[select],
                    ),
                    color=color,
                    label=label,
                    **kwargs
                )
            if ax_phase:
                angle_h = np.angle(h, deg=True)
                ax_phase.errorbar(
                    f[select],
                    angle_h[select],
                    (
                        angle_h[select]
                        - np.maximum((angle_h - 180 / np.pi / W)[select], -180),
                        np.minimum((angle_h + 180 / np.pi / W)[select], 180)
                        - angle_h[select],
                    ),
                    color=color,
                    label=label,
                    **kwargs
                )
                # ax_phase.fill_between(
                #    f[select],
                #    np.minimum((np.angle(h, deg = True) + 180/np.pi / W)[select], 180),
                #    np.maximum((np.angle(h, deg = True) - 180/np.pi / W)[select], -180),
                #    color = color,
                #    alpha = .2,
                #    **kwargs
                # )
        else:
            if ax_mag:
                db.mline = ax_mag.loglog(
                    fitter.F_Hz, abs(fitter.data), label=label, **kwargs
                )[0]

                if color is None:
                    color = db.mline.get_color()

            if ax_phase:
                _f, _h = angle_cutter(fitter.F_Hz, fitter.data, deg=True)
                db.pline = ax_phase.semilogx(
                    _f, _h, label=label, color=color, **kwargs
                )[0]

            db.min_data_err_lim = np.min(abs(h)) / 1.2
            db.max_data_err_lim = np.max(abs(h)) * 1.2

        midpt = (np.max(fitter.F_Hz) - np.min(fitter.F_Hz)) / 2
        if np.count_nonzero(fitter.F_Hz < midpt) < 0.65 * len(fitter.F_Hz):
            db.scale_log = False
        else:
            db.scale_log = True

        db.maxF_Hz = max(f)
        db.minF_Hz = min(f)
        select = f > 0
        db.minF_HzNZ = min(f[select])

        data_list.append(db)
        return db

    return plot


def plot_augment_fit(
    self,
    ax_mag,
    ax_phase,
    fit_list,
    with_error=False,
    F_Hz=None,
    label=args.UNSPEC,
):
    label_default = args.argscan(
        label, self.label_fit, "fit ({nP}P, {nZ}Z, Rsq ={Rsq:.3e})"
    )

    def plot(fitter, label=None, db_ref=None, color=None, **kwargs):
        if label is None:
            label = label_default
        label = label.format(
            nP=len(fitter.ZPKrep.poles),
            nZ=len(fitter.ZPKrep.zeros),
            Rsq=fitter.residuals_average,
        )

        if color is None and db_ref is not None:
            color = db_ref.mline.get_color()

        if F_Hz is not None:
            f = F_Hz
            h = fitter.xfer_eval(f)
        else:
            h = fitter.xfer_fit
            f = fitter.F_Hz
        db = Bunch()

        db.mline = ax_mag.loglog(f, abs(h), label=label, color=color, **kwargs)[0]

        if color is None:
            color = db.mline.get_color()

        if ax_phase:

            _f, _h = angle_cutter(f, h, deg=True)
            db.pline = ax_phase.semilogx(_f, _h, label=label, color=color, **kwargs)

        if with_error and isinstance(fitter, fitters_ZPK.MultiReprFilterZ):
            select = f > 0
            ax_mag.fill_between(
                f[select],
                np.maximum(
                    (abs(h) / (1 + fitter.xfer_fit_error.mag_rel2))[select],
                    np.min(abs(h)) / 10,
                ),
                np.minimum(
                    (abs(h) * (1 + fitter.xfer_fit_error.mag_rel2))[select],
                    np.max(abs(h)) * 10,
                ),
                color=db.mline.get_color(),
                alpha=0.3,
            )
            if ax_phase:
                ax_phase.fill_between(
                    f[select],
                    np.minimum(
                        (
                            np.angle(h, deg=True)
                            + 180 / np.pi * fitter.xfer_fit_error.phase
                        )[select],
                        180,
                    ),
                    np.maximum(
                        (
                            np.angle(h, deg=True)
                            - 180 / np.pi * fitter.xfer_fit_error.phase
                        )[select],
                        -180,
                    ),
                    color=db.pline[0].get_color(),
                    alpha=0.3,
                )
                db.min_fit_err_lim = np.min(abs(h)) * max(
                    1 / 3, np.min(1 / (1 + fitter.xfer_fit_error.mag_rel2)[select])
                )
                db.max_fit_err_lim = np.max(abs(h)) * min(
                    3, np.min((1 + fitter.xfer_fit_error.mag_rel2)[select])
                )
        else:
            db.min_fit_err_lim = np.min(abs(h)) / 1.2
            db.max_fit_err_lim = np.max(abs(h)) * 1.2
        db.median_fit_lim = np.median(abs(h))

        midpt = (np.max(fitter.F_Hz) - np.min(fitter.F_Hz)) / 2
        if np.count_nonzero(fitter.F_Hz < midpt) < 0.65 * len(fitter.F_Hz):
            db.scale_log = False
        else:
            db.scale_log = True

        db.maxF_Hz = max(f)
        db.minF_Hz = min(f)
        select = f > 0
        db.minF_HzNZ = min(f[select])

        fit_list.append(db)
        return db

    return plot


def plot_augment_residuals(
    self,
    ax_mag,
    ax_phase,
    fit_list,
    label=args.UNSPEC,
):
    label_default = args.argscan(
        label, self.label_fit, "fit ({nP}P, {nZ}Z, Rsq ={Rsq:.3e})"
    )

    def plot(
        fitter, label=None, db_ref=None, color=None, markersize=3, marker=".", **kwargs
    ):
        if label is None:
            label = label_default
        label = label.format(
            nP=len(fitter.ZPKrep.poles),
            nZ=len(fitter.ZPKrep.zeros),
            Rsq=fitter.residuals_average,
        )

        if color is None and db_ref is not None:
            color = db_ref.mline.get_color()

        r = fitter.residuals_preferred
        f = fitter.F_Hz
        db = Bunch()

        if markersize is not None:
            kwargs["s"] = markersize
        db.mline = ax_mag.scatter(
            f, abs(r), label=label, color=color, marker=marker, **kwargs
        )
        ax_mag.set_yscale("log_zoom")

        if color is None:
            color = db.mline.get_color()

        if ax_phase:
            _f, _r = angle_cutter(f, r, deg=True)
            db.pline = ax_phase.scatter(
                _f, _r, label=label, color=color, marker=marker, **kwargs
            )

        midpt = (np.max(f) - np.min(f)) / 2
        if np.count_nonzero(f < midpt) < 0.65 * len(f):
            db.scale_log = False
        else:
            db.scale_log = True

        db.maxF_Hz = max(f)
        db.minF_Hz = min(f)
        select = f > 0
        db.minF_HzNZ = min(f[select])

        db.min_res_err_lim = np.min(abs(r)) / 1.2
        db.max_res_err_lim = np.max(abs(r)) * 1.2
        db.median_res_lim = np.median(abs(r))
        fit_list.append(db)
        return db

    return plot


def plot_augment_zp(self, ax):
    def plot(fitter):
        BWs = []
        ZPKrep = ZPKrep2Sf(fitter.ZPKrep)
        p = ZPKrep.poles.fullplane
        F = p.imag
        BW = p.real
        select = BW < 0
        ax.errorbar(
            F[select],
            -BW[select],
            xerr=-BW[select],
            color=self.poleID_color,
            lw=0.5,
            alpha=0.5,
            marker="None",
            ls="None",
        )
        ax.scatter(
            F[select],
            -BW[select],
            s=10,
            facecolor=self.poleID_color,
            linewidths=0.5,
            # edgecolors = 'black',
        )
        BWs.extend(-BW[select])
        ax.errorbar(
            F[~select],
            BW[~select],
            xerr=BW[~select],
            color=self.poleOD_color,
            lw=0.5,
            alpha=0.5,
            marker="None",
            ls="None",
        )
        ax.scatter(
            F[~select],
            BW[~select],
            s=10,
            facecolor=self.poleOD_color,
            linewidths=0.5,
            # edgecolors = 'black',
        )
        BWs.extend(BW[~select])

        z = ZPKrep.zeros.fullplane
        F = z.imag
        BW = z.real
        select = BW < 0
        ax.errorbar(
            F[select],
            -BW[select],
            xerr=-BW[select],
            color=self.zeroID_color,
            lw=0.5,
            alpha=0.5,
            marker="None",
            ls="None",
        )
        ax.scatter(
            F[select],
            -BW[select],
            s=10,
            facecolor=self.zeroID_color,
            linewidths=0.5,
            # edgecolors = 'black',
            marker="x",
        )
        BWs.extend(-BW[select])
        ax.errorbar(
            F[~select],
            BW[~select],
            xerr=BW[~select],
            color=self.zeroOD_color,
            lw=0.5,
            alpha=0.5,
            marker="None",
            ls="None",
        )
        ax.scatter(
            F[~select],
            BW[~select],
            s=10,
            facecolor=self.zeroOD_color,
            linewidths=0.5,
            # edgecolors = 'black',
            marker="x",
        )
        BWs.extend(BW[~select])
        BWs = np.asarray(BWs)

        if len(z) > 0 or len(p) > 0:
            select = BWs > 0
            BWs = BWs[select]
            ax.set_ylim(np.max(BWs), np.min(BWs))
            ax.set_yscale("log_zoom")
        return

    return plot
