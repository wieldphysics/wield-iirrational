# -*- coding: utf-8 -*-
"""
"""


import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt

from ..utilities.mpl import generate_ax
from . import utilities
from .utilities import ZPKrep2Sf

def plot_ZP(
    self,
    fitter,
    db_nyquist = True,
    gs_base    = gridspec.GridSpec(1, 1)[0],
    fig        = None,
    horizontal = True,
    wspacing   = 0.25,
    hspace     = 0.2,
    ax_bunch   = None,
):
    axB = generate_ax(ax_bunch)

    if fig is None:
        fig = plt.figure()
        if horizontal:
            fig.set_size_inches(12, 5)
        else:
            fig.set_size_inches(5, 12)
    axB.fig = fig

    if horizontal:
        N_horizontal = 3
        N_vertical = 1
    else:
        N_horizontal = 1
        N_vertical = 3

    gs_DC = gridspec.GridSpecFromSubplotSpec(
        N_vertical, N_horizontal,
        subplot_spec = gs_base,
        wspace = wspacing, hspace = hspace,
        #height_ratios = height_ratio_list,
        #width_ratios = width_ratios,
    )

    if horizontal:
        axB.ax_poles = axB.fig.add_subplot(gs_DC[0, 0], polar = True)
        axB.ax_dual = axB.fig.add_subplot(gs_DC[0, 1], polar = True)
        axB.ax_zeros = axB.fig.add_subplot(gs_DC[0, 2], polar = True)
    else:
        axB.ax_poles = axB.fig.add_subplot(gs_DC[0, 0], polar = True)
        axB.ax_dual = axB.fig.add_subplot(gs_DC[1, 0], polar = True)
        axB.ax_zeros = axB.fig.add_subplot(gs_DC[2, 0], polar = True)
    axB.ax_poles.set_aspect('equal', 'box', anchor = 'C')
    axB.ax_dual.set_aspect('equal', 'box', anchor = 'C')
    axB.ax_zeros.set_aspect('equal', 'box', anchor = 'C')
    axB.ax_poles.set_rlabel_position(225)
    axB.ax_dual.set_rlabel_position(225)
    axB.ax_zeros.set_rlabel_position(225)

    def plot_poles(ax):
        p = fitter.poles_full
        r = abs(p)
        theta = np.angle(p)
        select = abs(p) < 1
        if db_nyquist:
            F_largest_Hz = np.max(fitter.F_Hz)
            r_pos = np.maximum(-10*np.log10(utilities.BWz(r[select], fitter.F_nyquist_Hz / F_largest_Hz)), 0)
            r_neg = np.maximum(-10*np.log10(utilities.BWz_neg(r[~select], fitter.F_nyquist_Hz / F_largest_Hz)), 0)
        else:
            r_pos = r[select]
            r_neg = r[~select]

        ax.scatter(
            theta[select],
            r_pos,
            s = 15,
            facecolor = self.poleID_color,
            linewidths = .3,
            edgecolors = 'black',
            label = 'stable(<1)',
        )
        ax.scatter(
            theta[~select],
            r_neg,
            s = 15,
            facecolor = self.poleOD_color,
            linewidths = .3,
            edgecolors = 'black',
            label = 'unstable(>1)',
        )
        return max(np.concatenate([[-10], r_pos, r_neg]))

    def plot_zeros(ax):
        z = fitter.zeros_full
        r = abs(z)
        theta = np.angle(z)
        select = abs(z) < 1
        if db_nyquist:
            F_largest_Hz = np.max(fitter.F_Hz)
            r_pos = np.maximum(-10*np.log10(utilities.BWz(r[select], fitter.F_nyquist_Hz / F_largest_Hz)), 0)
            r_neg = np.maximum(-10*np.log10(utilities.BWz_neg(r[~select], fitter.F_nyquist_Hz / F_largest_Hz)), 0)
        else:
            r_pos = r[select]
            r_neg = r[~select]

        ax.scatter(
            theta[select],
            r_pos,
            s = 15,
            facecolor = self.zeroID_color,
            linewidths = .5,
            #edgecolors = 'black',
            marker = 'x',
            label = 'mindelay (<1)',
        )
        ax.scatter(
            theta[~select],
            r_neg,
            s = 15,
            facecolor = self.zeroOD_color,
            linewidths = .5,
            #edgecolors = 'black',
            marker = 'x',
            label = 'delayed (>1)',
        )
        return max(np.concatenate([r_pos, r_neg]))

    r_max = plot_poles(axB.ax_poles)
    plot_poles(axB.ax_dual)
    r_max = max(plot_zeros(axB.ax_dual), r_max)
    plot_zeros(axB.ax_zeros)
    axB.ax_poles.set_title('poles')
    axB.ax_dual.set_title('dual')
    axB.ax_zeros.set_title('zeros')
    axB.ax_poles.set_ylim(0, (r_max // 3 + 1) * 3)
    axB.ax_dual.set_ylim(0, (r_max // 3 + 1) * 3)
    axB.ax_zeros.set_ylim(0, (r_max // 3 + 1) * 3)
    #ax1.axhline(40, color = 'purple')
    #axB.ax_dual.axhline(40, color = 'purple')
    #axB.ax_zeros.axhline(40, color = 'purple')

    num = 9
    axB.ax_poles.set_xticks(np.linspace(0, np.pi, num))
    axB.ax_poles.set_xticklabels(["{0:.1f}".format(f) for f in np.linspace(0, fitter.F_nyquist_Hz, num)])
    axB.ax_dual.set_xticks(np.linspace(0, np.pi, num))
    axB.ax_dual.set_xticklabels(["{0:.1f}".format(f) for f in np.linspace(0, fitter.F_nyquist_Hz, num)])
    axB.ax_zeros.set_xticks(np.linspace(0, np.pi, num))
    axB.ax_zeros.set_xticklabels(["{0:.1f}".format(f) for f in np.linspace(0, fitter.F_nyquist_Hz, num)])
    #c = patches.Circle((0, 0), 1, color = 'black', fill = False)
    #ax.add_artist(c)
    axB.ax_poles.legend(loc = 'lower right', fontsize = 8)
    axB.ax_zeros.legend(loc = 'lower right', fontsize = 8)
    axB.ax_poles.grid(b=True)
    axB.ax_dual.grid(b=True)
    axB.ax_zeros.grid(b=True)

    axB.finalize()
    return axB

def plot_ZP_S(
    self,
    fitter,
    db_nyquist = True,
    gs_base    = gridspec.GridSpec(1, 1)[0],
    fig        = None,
    horizontal = True,
    wspacing   = 0.25,
    hspace     = 0.2,
    ax_bunch   = None,
):
    axB = generate_ax(ax_bunch)

    if fig is None:
        fig = plt.figure()
        if horizontal:
            fig.set_size_inches(12, 5)
        else:
            fig.set_size_inches(5, 12)
    axB.fig = fig

    if horizontal:
        N_horizontal = 3
        N_vertical = 1
    else:
        N_horizontal = 1
        N_vertical = 3

    gs_DC = gridspec.GridSpecFromSubplotSpec(
        N_vertical, N_horizontal,
        subplot_spec = gs_base,
        wspace = wspacing, hspace = hspace,
        #height_ratios = height_ratio_list,
        #width_ratios = width_ratios,
    )

    if horizontal:
        axB.ax_poles = axB.fig.add_subplot(gs_DC[0, 0])
        axB.ax_dual = axB.fig.add_subplot(gs_DC[0, 1])
        axB.ax_zeros = axB.fig.add_subplot(gs_DC[0, 2])
    else:
        axB.ax_poles = axB.fig.add_subplot(gs_DC[0, 0])
        axB.ax_dual = axB.fig.add_subplot(gs_DC[1, 0])
        axB.ax_zeros = axB.fig.add_subplot(gs_DC[2, 0])
    axB.ax_poles.yaxis.tick_right()
    axB.ax_dual.yaxis.tick_right()
    axB.ax_zeros.yaxis.tick_right()

    axB.ax_poles.yaxis.set_label_position("right")
    axB.ax_dual.yaxis.set_label_position("right")
    axB.ax_zeros.yaxis.set_label_position("right")

    def plot_poles(ax, ZPKrep):
        p = ZPKrep.poles.fullplane
        F_Hz = p.imag
        BW_Hz = p.real
        select = p.real < 0

        ax.scatter(
            -BW_Hz[select],
            F_Hz[select],
            s = 15,
            facecolor = self.poleID_color,
            linewidths = .3,
            edgecolors = 'black',
            label = 'stable(re<0)',
        )
        ax.scatter(
            BW_Hz[~select],
            F_Hz[~select],
            s = 15,
            facecolor = self.poleOD_color,
            linewidths = .3,
            edgecolors = 'black',
            label = 'unstable(re>0)',
        )
        return BW_Hz

    def plot_zeros(ax, ZPKrep):
        z = ZPKrep.zeros.fullplane
        F_Hz = z.imag
        BW_Hz = z.real
        select = z.real < 0

        ax.scatter(
            -BW_Hz[select],
            F_Hz[select],
            s = 15,
            facecolor = self.zeroID_color,
            linewidths = .5,
            #edgecolors = 'black',
            marker = 'x',
            label = 'mindelay (re<0)',
        )
        ax.scatter(
            BW_Hz[~select],
            F_Hz[~select],
            s = 15,
            facecolor = self.zeroOD_color,
            linewidths = .5,
            #edgecolors = 'black',
            marker = 'x',
            label = 'delayed (re>0)',
        )
        return BW_Hz

    def x_limits(*r_vals):
        r_vals = np.concatenate(r_vals)
        if len(r_vals) == 0:
            return 1, .1
        rmin = min(abs(r_vals))
        rmax = max(abs(r_vals))
        return rmax * 3, rmin / 3

    ZPKrep = ZPKrep2Sf(fitter.ZPKrep)
    r_p = plot_poles(axB.ax_poles, ZPKrep)
    plot_poles(axB.ax_dual, ZPKrep)
    r_z = plot_zeros(axB.ax_dual, ZPKrep)
    plot_zeros(axB.ax_zeros, ZPKrep)
    axB.ax_poles.set_ylabel('poles Freq [Hz]')
    axB.ax_dual.set_ylabel('dual Freq [Hz]')
    axB.ax_zeros.set_ylabel('zeros Freq [Hz]')

    if len(r_p) > 0:
        axB.ax_poles.set_xscale('log_zoom')

    axB.ax_poles.set_xlim(x_limits(r_p))
    if len(r_p) > 0 or len(r_z) > 0:
        axB.ax_dual.set_xscale('log_zoom')
    axB.ax_dual.set_xlim(x_limits(r_p, r_z))
    if len(r_z) > 0:
        axB.ax_zeros.set_xscale('log_zoom')
    axB.ax_zeros.set_xlim(x_limits(r_z))

    axB.ax_poles.legend(loc = 'lower left', fontsize = 8)
    axB.ax_zeros.legend(loc = 'lower left', fontsize = 8)
    axB.ax_poles.grid(b=True)
    axB.ax_dual.grid(b=True)
    axB.ax_zeros.grid(b=True)

    axB.finalize()
    return axB
