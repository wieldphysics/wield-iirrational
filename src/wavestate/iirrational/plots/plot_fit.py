# -*- coding: utf-8 -*-
"""
"""


import numpy as np
from matplotlib import gridspec

from ..utilities.np import logspaced
from ..utilities.mpl import generate_stacked_plot_ax
from ..utilities import args

from .utilities import ZPKrep2Sf
from . import plot_augment


def plot_fit(
        self,
        fitter,
        plot_zp         = True,
        plot_past_data  = False,
        plot_to_nyquist = True,
        plot_data_scale     = 10,
        plot_data_scale_low = 2,
        with_error      = True,
        with_data_error = True,
        plot_data       = True,
        N_points        = 150,
        gs_base         = gridspec.GridSpec(1, 1)[0],
        fig             = None,
        ax_bunch        = None,
        xscale          = None,
        plot_now        = True,
        label_data      = args.UNSPEC,
        label_fit       = args.UNSPEC,
):
    use_phase = True
    label_data = args.argscan(label_data, self.label_data, 'data')
    label_fit  = args.argscan(label_fit, self.label_fit, 'fit ({nP}P, {nZ}Z, Rsq ={Rsq:.2e})')
    #if fitter.magsq_data:
    #    use_phase = False
    if fitter.phase_missing:
        use_phase = False

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

    ZPKrepSf = ZPKrep2Sf(fitter.ZPKrep)
    if len(ZPKrepSf.poles) + len(ZPKrepSf.zeros) == 0:
        plot_zp = False
    if plot_zp:
        axB.mag2 = axB.mag.twinx()

    data_list = []
    if plot_data and fitter.data is not None:
        axB.plot_data = plot_augment.plot_augment_data(
            self,
            ax_mag     = axB.mag,
            ax_phase   = axB.phase,
            data_list  = data_list,
            with_error = with_error,
            label      = label_data,
        )
        if plot_now:
            axB.plot_data(fitter)

    fit_list = []
    if plot_past_data:
        intersperse = min((len(fitter.F_Hz)-1) // N_points + 1, 1)
        if fitter.F_nyquist_Hz is not None:
            if plot_to_nyquist:
                f = np.sort(np.concatenate([
                    np.linspace(0, fitter.F_nyquist_Hz, N_points//2),
                    logspaced(np.min(fitter.F_Hz[fitter.F_Hz > 0]) / plot_data_scale_low, fitter.F_nyquist_Hz, N_points//2),
                    fitter.F_Hz[intersperse//2::intersperse],
                ]))
            else:
                F_max_Hz = np.max(fitter.F_Hz) * plot_data_scale
                f = np.sort(np.concatenate([
                    np.linspace(0, F_max_Hz, N_points//2),
                    logspaced(np.min(fitter.F_Hz[fitter.F_Hz > 0]) / plot_data_scale_low, F_max_Hz, N_points//2),
                    fitter.F_Hz[intersperse//2::intersperse],
                ]))
        else:
            F_max_Hz = np.max(fitter.F_Hz) * plot_data_scale
            f = np.sort(np.concatenate([
                np.linspace(0, F_max_Hz, N_points//2),
                logspaced(np.min(fitter.F_Hz[fitter.F_Hz > 0]) / plot_data_scale_low, F_max_Hz, N_points//2),
                fitter.F_Hz[intersperse//2::intersperse],
            ]))
        axB.plot_fit = plot_augment.plot_augment_fit(
            self,
            ax_mag     = axB.mag,
            ax_phase   = axB.phase,
            with_error = with_error,
            label      = label_fit,
            F_Hz       = f,
            fit_list   = fit_list,
        )
    else:
        axB.plot_fit = plot_augment.plot_augment_fit(
            self,
            ax_mag     = axB.mag,
            ax_phase   = axB.phase,
            with_error = with_error,
            label      = label_fit,
            fit_list   = fit_list,
        )
    if plot_now:
        axB.plot_fit(fitter)

    if plot_zp:
        axB.plot_zp = plot_augment.plot_augment_zp(self, ax = axB.mag2)
        if plot_now:
            axB.plot_zp(fitter)
    else:
        axB.plot_zp = lambda fitter : None

    axB.mag.set_ylabel('Magnitude')
    if axB.phase:
        axB.phase.set_ylabel('Phase [deg]')
    axB.ax_bottom.set_xlabel('Frequency [Hz]')

    def finalize():
        min_lims = []
        max_lims = []
        for fitB in data_list:
            min_err_lim = min(
                fitB.min_data_err_lim,
                max(fitB.min_data_err_lim / 30, fitB.min_data_err_lim * 1)
            )
            min_lims.append(min_err_lim)
            max_err_lim = max(
                fitB.max_data_err_lim,
                min(fitB.max_data_err_lim * 30, fitB.max_data_err_lim / 1)
            )
            max_lims.append(max_err_lim)
        for fitB in fit_list:
            min_err_lim = min(
                fitB.median_fit_lim,
                fitB.min_fit_err_lim,
                max(fitB.min_fit_err_lim / 30, fitB.min_fit_err_lim * 1)
            )
            min_lims.append(min_err_lim)
            max_err_lim = max(
                fitB.median_fit_lim,
                fitB.max_fit_err_lim,
                min(fitB.max_fit_err_lim * 30, fitB.max_fit_err_lim / 1)
            )
            max_lims.append(max_err_lim)
        min_lims = np.asarray(min_lims)
        max_lims = np.asarray(max_lims)
        axB.mag.set_ylim(
            np.min(min_lims[np.isfinite(min_lims)]),
            np.max(max_lims[np.isfinite(max_lims)])
        )

        scales = []
        mins   = []
        minsNZ = []
        maxs   = []
        for fitB in fit_list + data_list:
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
            ncol = 1,
            fontsize = 8,
            loc = 'best'
        )
    axB.finalizers.append(finalize)
    return axB
