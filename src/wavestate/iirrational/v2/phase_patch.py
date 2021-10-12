# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
import declarative

from ..utilities import ensure_aid

from ..utilities import np as np_utils
from .. import representations
from .. import fitters_ZPK

from . import algorithms

def fitter_log_mag_generate(aid):
    #TODO
    return

def fitter_stabilized_log_mag_generate(aid):
    #TODO
    return

def root_stabilize(aid):
    z_c = np.copy(aid.fitter.zeros.c)
    z_r = np.copy(aid.fitter.zeros.r)
    z_c.real[z_c.real > 0] *= -1
    z_r.real[z_r.real > 0] *= -1

    p_c = np.copy(aid.fitter.poles.c)
    p_r = np.copy(aid.fitter.poles.r)
    p_c.real[p_c.real > 0] *= -1
    p_r.real[p_r.real > 0] *= -1

    aid.fitter_update(
        aid.fitter.regenerate(
            coding_map = fitters_ZPK.codings_s.coding_maps.SOSsafe,
            z_c = z_c,
            z_r = z_r,
            p_c = p_c,
            p_r = p_r,
        ),
        representative = False,
    )
    return

def phase_patch(aid):
    fitter_mag = aid.fitter.copy()
    fitter_mag.residuals_log_im_scale = 0
    with fitter_mag.with_codings_only([fitter_mag.gain_coding]):
        fitter_mag.optimize_NM(residuals_type = 'log')
    fitter_mag.optimize(
        aid            = aid,
        residuals_type = 'log',
    )
    res_avg = fitter_mag.residuals_average
    aid.log_debug(5, 'log residuals average', res_avg)

    Wsort = np.sort(aid.fitter.W)
    Wthresh = Wsort[int(len(aid.fitter.W) * .2)]
    #TODO decide about this
    Wselect = aid.fitter.W >= 0 * Wthresh / 2
    data = (fitter_mag.data / fitter_mag.xfer_fit)
    #data = data / abs(data)
    xres = representations.ZPKwData(
        data = data[Wselect],
        W = aid.fitter.W[Wselect],
        F_Hz = aid.fitter.F_Hz[Wselect],
        F_nyquist_Hz = None,
    )

    fitter_X = fitters_ZPK.ZPKrep2MRF(
        ZPKrep = xres,
        delay_s_min = 0,
        delay_s_max = 0,
    )
    fitter_X.gain = 1
    #only use the first fraction of the data, the latter data is the most influenced
    #by delay and by bad fitting
    algorithms.sign_check_flip(fitter_X, max_N = len(fitter_X.F_Hz) // 3 + 3)
    fitter_X.optimize(residuals_type = 'zeros')
    fitter_X.delay_coding = fitters_ZPK.codings_cmn.CodingDelay(fitter_X)
    fitter_X.optimize(residuals_type = 'log')

    if aid.hint('trust_SNR'):
        weight_scale = 1
    else:
        with aid.log_heading('sample variance (from magnitude)'):
            R = (fitter_X.data / fitter_X.xfer_fit)
            W = fitter_X.W
            meanWeight = (np.sum(W**4) / np.sum(W**2))**.5
            cRvar = (np.sum((abs(R) - 1)**2 * W**2) / np.sum(W**2))**.5
            aid.log_debug(5, 'variance', 1/cRvar, meanWeight)
            weight_scale_true = meanWeight * cRvar
            weight_scale = weight_scale_true**.5
            aid.log_alert(3, "Weight Scaling determined: ", weight_scale_true)
            aid.log_alert(3, "Weight Scaling Used ", weight_scale)

    phase_fixes = []
    with aid.log_heading('Phase patching'):
        shift_saves = []
        separations = list(np_utils.logspaced(
            3,
            len(fitter_X.F_Hz) / 2,
            int(1 + np.log(len(fitter_X.F_Hz)))
        ))
        check_second_pass = False
        while separations:
            sep = separations.pop()
            sep = int(sep)

            f = fitter_X.F_Hz
            f_d = (f[2:] - f[:-2]) / 2
            W = fitter_X.W[1:-1]
            R = (fitter_X.data / fitter_X.xfer_fit)[1:-1]
            f = f[1:-1]

            cR = np.cumsum(R/abs(R) * W**2 * f_d)
            cW = np.cumsum(W**2 * f_d)
            cFd = np.cumsum(f_d)

            shift_save = wavestate.bunch.Bunch()
            #distance_10
            shift = (cR[sep:] - cR[0:-sep]) / (cW[sep:] - cW[0:-sep])
            shift_err = weight_scale * ((cFd[sep:] - cFd[0:-sep]) / (cW[sep:] - cW[0:-sep])  / (2*sep))**.5
            f_local = (f[sep:] + f[:-sep])/2
            bw_count = 0
            for sec_L, sec_R in np_utils.generate_sections(shift + shift_err < 0 ):
                loc_Hz = (f_local[sec_R] + f_local[sec_L])/2
                BW_Hz = (f_local[sec_R] - f_local[sec_L]) / 2
                aid.log(6, "LOCATION ", loc_Hz)
                aid.log(6, "BW: ", BW_Hz)

                try:
                    select2BW = (f_local < loc_Hz + 2*BW_Hz) & (f_local > loc_Hz - 2*BW_Hz)
                    Imax = np.max(shift.imag[select2BW])
                    Imin = np.min(shift.imag[select2BW])
                    Isep = (Imax - Imin)
                    Icmn = (Imax + Imin)/2
                    aid.log(8, "Isep", Isep)
                    aid.log(8, "Icmn", Icmn)

                    if abs(Isep) + shift_err[sec_R] + shift_err[sec_L] < 0.75:
                        aid.log_debug(8, "Skipping Pair")
                        continue
                except ValueError:
                    aid.log_debug(8, "Skipping Pair, phase patch missing")
                    continue

                bw_count += 1
                if check_second_pass:
                    continue

                coding = fitters_ZPK.codings_cmn.CodingDelayPair(fitter_X)
                coding.F_Hz = loc_Hz
                coding.BW_Hz = BW_Hz
                fitter_X.num_codings.append(coding)
                fitter_X.codings_revision += 1

            if not check_second_pass:
                shift_save.sep = sep
                shift_save.f_local = f_local
                shift_save.shift = shift
                shift_save.shift_err = shift_err
                shift_saves.append(shift_save)
                fitter_X.optimize()
                #switch to second pass mode check
                check_second_pass = True
                #add the separation back for the second pass
                separations.append(sep)
            else:
                if bw_count > 0:
                    aid.log_warn(
                        1,
                        """
                        Phase patching did not remove all unstable roots in
                        a single pass. It is not programmed to accommodate this
                        failure mode. Unstable zeros (or poles) may be wrong.
                        """)
                check_second_pass = False

            phase_fix_about = wavestate.bunch.Bunch()
            phase_fix_about.shift_saves = shift_saves
            phase_fix_about.fitter_X = fitter_X.copy()
            phase_fixes.append(phase_fix_about)

    def investigate_phase_patch_plot(
            self,
            fname = None,
            xscale = 'log_zoom',
            **kwargs
    ):
        """
        Generates a plot of the mag normalized fit. Returns a wavestate.bunch.Bunch object which
        stores the axes and figure of the plot.
        """
        from .. import plots
        axB = plots.plot_fitter_flag(
            fitter = phase_fixes[0].fitter_X,
            xscale = xscale,
            **kwargs
        )
        for pf in phase_fixes[1:]:
            axB.plot_fit(pf.fitter_X)
        if fname is not None:
            axB.save(fname)
        return axB
    aid.investigations['investigate_phase_patch_plot'] = investigate_phase_patch_plot

    def investigate_phase_patch_stats_plot(
            resaid,
            fname = None,
            xscale = 'log_zoom'
    ):
        """
        Boy does this plot need documentation.
        """
        #TODO, put into a function for the annotater
        from IIRrational.utilities.mpl import mplfigB
        axB = mplfigB(Ncols = 2, Nrows=2)
        for shift in shift_saves:
            axB.ax0_0.plot(shift.f_local, shift.shift.real)
            axB.ax1_0.plot(shift.f_local, shift.shift.imag)
            axB.ax0_1.plot(shift.f_local, shift.shift_err)
        axB.ax0_0.set_xscale(xscale)
        axB.ax1_0.set_xscale(xscale)
        axB.ax0_1.set_xscale(xscale)
        axB.finalize()
        if fname is not None:
            axB.save(fname)
        return axB
    aid.investigations['investigate_phase_patch_stats_plot'] = investigate_phase_patch_stats_plot

    if len(fitter_X.num_codings) > 10:
        aid.log_warn(
            3, ("""
            Phase patching found an excessive number of phase loops. Rational fitting probably didn't perform well... Skipping phase patching.
            Check --mode=rational to see how bad that fit was. Use --poles and --zeros suggestions to improve it.
            """)
        )
        return aid.fitter

    extra_p_c = []
    extra_z_c = []
    for coding in fitter_X.num_codings:
        i = coding.F_Hz
        r = coding.BW_Hz
        if r < 0:
            if aid.hint('never_unstable_poles'):
                aid.log_warn(
                    3, ("""
                    High Confidence that an unstable pole exists at {}+{}i Hz!
                    Not adding it to the filter, due to 'never_unstable_poles'
                    """).format(-r, i)
                )
                continue
            else:
                aid.log_warn(
                    3, ("""
                    High Confidence that an unstable pole exists at {}+{}i Hz!
                    Adding it to the filter. Prevent this with 'never_unstable_poles'
                    """).format(-r, i)
                )
        else:
            if aid.hint('never_unstable_zeros'):
                aid.log_warn(
                    8, ("""
                    High Confidence that an unstable zero exists at {}+{}i Hz
                    Not adding it to the filter, due to 'never_unstable_zeros'
                    """).format(+r, i)
                )
                continue
            else:
                aid.log_warn(
                    5, ("""
                    High Confidence that an unstable zero exists at {}+{}i Hz
                    Adding it to the filter. Prevent this with 'never_unstable_zeros'
                    """).format(+r, i)
                )
        extra_p_c.append(-r + 1j*i)
        extra_z_c.append(+r + 1j*i)

    extra_p_r = []
    extra_z_r = []
    ord_diff = -int(4 * fitter_X.delay_s * fitter_X.F_max_Hz)
    reldeg_max = aid.hint('relative_degree_max')
    reldeg_min = aid.hint('relative_degree_min')
    old_order = aid.fitter_orders().reldeg
    new_order = old_order + ord_diff
    if ord_diff < 0:
        if reldeg_min is not None and new_order < reldeg_min:
            ord_diff = old_order - reldeg_min
        extra_p_r = [-aid.fitter.F_max_Hz] * -ord_diff
    else:
        if reldeg_max is not None and new_order > reldeg_max:
            ord_diff = reldeg_max - old_order
        extra_z_r = [-aid.fitter.F_max_Hz] * ord_diff
    gain_adjust = aid.fitter.F_max_Hz**ord_diff

    #TODO add restrictions on the order increase-decrease
    fitter_new = aid.fitter.regenerate(
        z_c = list(fitter_mag.zeros.c) + extra_z_c,
        p_c = list(fitter_mag.poles.c) + extra_p_c,
        z_r = list(fitter_mag.zeros.r) + extra_z_r,
        p_r = list(fitter_mag.poles.r) + extra_p_r,
    )
    fitter_new.gain = fitter_mag.gain * gain_adjust
    #TODO make this set to the specified delay
    fitter_new.delay_s = 0

    with fitter_new.with_codings_only([fitter_new.gain_coding]):
        fitter_new.optimize()
    algorithms.optimize_anneal(aid, fitter_new)
    aid.fitter_update(
        fitter_new,
        representative = True,
    )
    return fitter_new

