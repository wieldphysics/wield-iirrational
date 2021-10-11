# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
from matplotlib import gridspec

from ..utilities.np import logspaced
from ..utilities.mpl import generate_stacked_plot_ax

def plot_rel_comparison(
    self,
    fitter,
    fitter_ref,
    N_points        = 1000,
    gs_base         = gridspec.GridSpec(1, 1)[0],
    fig             = None,
    ax_bunch        = None,
    xscale          = None,
):
    if ax_bunch is None:
        axB = generate_stacked_plot_ax(
            [
                ('rel', True),
            ],
            height_ratios = dict(
                mag = 1,
            ),
            heights_phys_in_default = 1,
            xscales = 'log',
            gs_base = gs_base,
            fig     = fig,
            ax_bunch = ax_bunch,
        )
    else:
        axB = ax_bunch

    axB.phase = axB.rel.twinx()
    f_min = np.max(fitter.F_Hz)
    f = np.sort(
        np.concatenate([
            fitter.F_Hz,
            np.linspace(f_min, fitter.F_nyquist_Hz, N_points//2),
            logspaced(f_min, fitter.F_nyquist_Hz, N_points//4),
        ]))
    f = f[f > 0]
    h = fitter.xfer_eval(f)
    h_ref = fitter_ref.xfer_eval(f)

    rel = h / h_ref

    mline = axB.rel.semilogx(
        f,
        20 * np.log10(abs(rel)),
        color = 'black',
        label = 'rel mag [db]'
    )
    pline = axB.phase.semilogx(
        f,
        np.angle(rel),
        label = 'fit phase [rad]',
        color = 'purple',
    )

    axB.rel.legend(
        handles = [mline[0], pline[0]],
        ncol = 2,
        fontsize = 8,
        loc = 'lower right'
    )
    axB.rel.set_xlabel('Frequency [Hz]')
    axB.rel.set_ylabel('DB difference')
    axB.phase.set_ylabel('Rad Difference')

    axB.finalize()
    return axB
