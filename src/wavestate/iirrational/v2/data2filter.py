# -*- coding: utf-8 -*-
"""
"""

import numpy as np

from ..utilities import np as util_np
from .. import fitters_ZPK
from .. import representations
from ..utilities import np as np_utils

from . import arguments
from .arguments import (
    grab_kwargs,
    grab_kwarg_hints,
)
from . import order_reduce
from . import algorithms
from . import phase_patch
from . import fit_aid
from . import results_aid_adv
from . import rational_fits
from . import SNR_adjustments
from . import autogen_docstr

def data2filter(*args, **kw):
    """
    """
    kwput = dict()
    kw, ZPKrep = arguments.kw_ZPKrep_build(args, kw)

    #pop all advanced hints
    hints = kw.pop('hints', {})
    aid = fit_aid.FitAid(hints = hints)

    #setup logging early so that error reporting preferences are done before any
    #errors are reported
    grab_kwarg_hints(aid, kw, arguments.logging.kw_hints, kwput = kwput)

    helparg = grab_kwargs(aid, kw, arguments.kw_hints, 'help')

    if helparg:
        raise NotImplementedError("Can't provide help functionality yet")

    #now pull the defaults from chain
    chain = grab_kwargs(aid, kw, arguments.kw_hints, 'chain')
    if chain is not None:
        raise NotImplementedError("WIP")

    grab_kwarg_hints(aid, kw, arguments.standardargs.kw_hints_pre, kwput = kwput)
    data          = grab_kwargs(aid, kw, arguments.kw_hints, 'data')
    F_Hz          = grab_kwargs(aid, kw, arguments.kw_hints, 'F_Hz')
    SNR           = grab_kwargs(aid, kw, arguments.kw_hints, 'SNR')
    COH           = grab_kwargs(aid, kw, arguments.kw_hints, 'coherence')
    SNR_phase_rel = grab_kwargs(aid, kw, arguments.kw_hints, 'SNR_phase_relative')
    emphasis      = grab_kwargs(aid, kw, arguments.kw_hints, 'emphasis')
    select        = grab_kwargs(aid, kw, arguments.kw_hints, 'select')
    F_nyquist_Hz  = grab_kwargs(aid, kw, arguments.kw_hints, 'F_nyquist_Hz', kwput = kwput)

    mode         = grab_kwargs(aid, kw, arguments.kw_hints, 'mode', kwput = kwput)

    ZPK   = grab_kwargs(aid, kw, arguments.kw_hints, 'ZPK', kwput = kwput)
    zeros = grab_kwargs(aid, kw, arguments.kw_hints, 'zeros', kwput = kwput)
    poles = grab_kwargs(aid, kw, arguments.kw_hints, 'poles', kwput = kwput)
    gain  = grab_kwargs(aid, kw, arguments.kw_hints, 'gain', kwput = kwput)

    if zeros is None and ZPK is not None:
        zeros = ZPK.zeros
    if zeros is None and ZPKrep is not None:
        zeros = ZPKrep.zeros
    if zeros is None:
        zeros = ()

    if poles is None and ZPK is not None:
        poles = ZPK.poles
    if poles is None and ZPKrep is not None:
        poles = ZPKrep.poles
    if poles is None:
        poles = ()

    if gain is None and ZPK is not None:
        gain = ZPK.gain
    if gain is None and ZPKrep is not None:
        gain = ZPKrep.gain
    if gain is None:
        gain = 1

    ZPK_overlay   = grab_kwargs(aid, kw, arguments.kw_hints, 'ZPK_overlay', kwput = kwput)
    zeros_overlay = grab_kwargs(aid, kw, arguments.kw_hints, 'zeros_overlay', kwput = kwput)
    poles_overlay = grab_kwargs(aid, kw, arguments.kw_hints, 'poles_overlay', kwput = kwput)

    if zeros_overlay is None and ZPK_overlay is not None:
        zeros_overlay = ZPK_overlay.zeros
    if zeros_overlay is None and ZPKrep is not None:
        zeros_overlay = ZPKrep.zeros_overlay
    if zeros_overlay is None:
        zeros_overlay = ()

    if poles_overlay is None and ZPK_overlay is not None:
        poles_overlay = ZPK_overlay.poles
    if poles_overlay is None and ZPKrep is not None:
        poles_overlay = ZPKrep.poles_overlay
    if poles_overlay is None:
        poles_overlay = ()

    if F_nyquist_Hz is not None:
        raise NotImplementedError("V2 currently only works in the Sf domain")

    if F_Hz is None and ZPKrep is not None:
        F_Hz = ZPKrep.F_Hz
    if data is None and ZPKrep is not None:
        data = ZPKrep.data
    if SNR is None and ZPKrep is not None:
        SNR = ZPKrep.W

    if SNR_phase_rel is None and ZPKrep is not None:
        SNR_phase_rel = ZPKrep.residuals_log_im_scale

    if SNR is not None:
        if isinstance(SNR, tuple):
            SNR = tuple(np_utils.broadcast_arrays_none(*SNR))
        else:
            SNR = (SNR,)

    if COH is not None:
        if isinstance(COH, tuple):
            COH = tuple(np_utils.broadcast_arrays_none(*COH))
        else:
            COH = (COH,)

        SNRs = []
        for coh in COH:
            SNRs.append((1/(1 - coh**0.5) - 1)**0.5)
        if SNR is None:
            SNR = tuple(SNRs)
        else:
            SNR = SNR + tuple(SNRs)

    #now combine all of the SNR arrays into a single one
    if SNR is not None:
        iSNR = 0
        for snr in SNR:
            iSNR = iSNR + snr**-2
        SNR = iSNR**-0.5

    if SNR_phase_rel is None:
        SNR_phase_rel = 1

    F_Hz, data, SNR, SNR_phase_rel, emphasis = np_utils.broadcast_arrays_none(
        F_Hz, data,
        SNR, SNR_phase_rel,
        emphasis
    )
    F_Hz, data, SNR, SNR_phase_rel, emphasis = np_utils.select_through_none(
        select, F_Hz, data,
        SNR, SNR_phase_rel,
        emphasis
    )

    SNR_min = grab_kwargs(aid, kw, arguments.adjustments.kw_hints, 'SNR_min', kwput = kwput)
    SNR_est_width = grab_kwargs(aid, kw, arguments.kw_hints, 'SNR_estimate_width', kwput = kwput)
    if SNR is None:
        if SNR_est_width is None or SNR_est_width == 0:
            aid.log_warn(3, (
                """
                SNR not specified, using unweighted "1". If data has many points
                with low SNR, the fit will be poor. Try 'SNR_est_width' to attempt
                an estimate using sample variance.
                """
            ).format())
            SNR = np.ones_like(F_Hz)
        else:
            aid.log_warn(3, (
                """
                Estimating SNR from sample variance with nearby points (SNR_est_width={} > 0).
                This technique works semi-OK, but could probably be much better..
                use the resulting fit to estimate the sample variance and generate
                improved SNR estimates, iterate.
                """
            ).format(SNR_est_width))
            from . import SNR_estimate
            SNR = SNR_estimate.SNR_estimate(F_Hz, data, width = SNR_est_width)
            if SNR_min is None:
                SNR_min = 1
            W_select = SNR > SNR_min
            if not np.all(W_select):
                F_Hz, data, SNR, SNR_phase_rel, emphasis = np_utils.select_through_none(W_select, F_Hz, data, SNR, SNR_phase_rel, emphasis)
                aid.log_warn(3, (
                    """{} SNR<{} element(s) dropped (of {}).
                    Too many low SNR elements confuses the rational nonparametric fitter.
                    """
                ).format(np.count_nonzero(~W_select), SNR_min, len(W_select)))
    else:
        if SNR_min is None:
            SNR_min = 0
        W_select = SNR > SNR_min
        if not np.all(W_select):
            F_Hz, data, SNR, SNR_phase_rel, emphasis = np_utils.select_through_none(W_select, F_Hz, data, SNR, SNR_phase_rel, emphasis)
            aid.log_warn(3, "{} SNR=0 element(s) dropped (of {})".format(np.count_nonzero(~W_select), len(W_select)))

    F_max_Hz = grab_kwargs(aid, kw, arguments.adjustments.kw_hints, 'F_max_Hz', kwput = kwput)
    if F_max_Hz is not None:
        select_Fmax = (F_Hz <= F_max_Hz)
        aid.log_warn(3, "{} element(s) dropped (of {}) due to F_max_Hz={}".format(np.count_nonzero(~select_Fmax), len(select_Fmax), F_max_Hz))
        F_Hz, data, SNR, SNR_phase_rel, emphasis = np_utils.select_through_none(select_Fmax, F_Hz, data, SNR, SNR_phase_rel, emphasis)

    F_min_Hz = grab_kwargs(aid, kw, arguments.adjustments.kw_hints, 'F_min_Hz', kwput = kwput)
    if F_min_Hz is not None:
        select_Fmax = (F_Hz >= F_min_Hz)
        aid.log_warn(3, "{} element(s) dropped (of {}) due to F_min_Hz={}".format(np.count_nonzero(~select_Fmax), len(select_Fmax), F_min_Hz))
        F_Hz, data, SNR, SNR_phase_rel, emphasis = np_utils.select_through_none(select_Fmax, F_Hz, data, SNR, SNR_phase_rel, emphasis)

    nan_select = np.isfinite(data)
    if not np.all(nan_select):
        F_Hz, data, SNR, SNR_phase_rel, emphasis = np_utils.select_through_none(nan_select, F_Hz, data, SNR, SNR_phase_rel, emphasis)
        aid.log_warn(3, "NaN's and/or inf's present in dataset. Dropping them.")


    ###################################### downsampling #######################
    downsampling_N = grab_kwargs(aid, kw, arguments.adjustments.kw_hints, 'downsample', kwput = kwput)
    downsampling_type = grab_kwargs(aid, kw, arguments.adjustments.kw_hints, 'downsample_type', kwput = kwput)
    downsampling_type = downsampling_type.lower()
    if downsampling_N is not None and downsampling_N < len(F_Hz):
        if downsampling_type == 'log':
            F_groups = util_np.logspaced(F_Hz[0], F_Hz[-1], downsampling_N)
        if downsampling_type in ['linear', 'lin']:
            F_groups = np.linspace(F_Hz[0], F_Hz[-1], downsampling_N)
        if downsampling_type in ['loglin', 'linlog']:
            F_groups = util_np.loglinspaced(F_Hz[0], F_Hz[-1], downsampling_N)

        idx_groups = np.searchsorted(F_Hz, F_groups)
        idx_pairs = list(zip(idx_groups[:-1], idx_groups[1:]))

        lemph = []
        lF_Hz  = []
        ldata = []
        lSNR = []
        lSNRpr = []
        for idx1, idx2 in idx_pairs:
            if idx1 == idx2:
                continue
            ldata.append(
                np.sum(data[idx1: idx2] * SNR[idx1: idx2]) / np.sum(SNR[idx1: idx2])
            )
            lF_Hz.append(
                np.sum(F_Hz[idx1: idx2] * SNR[idx1: idx2]) / np.sum(SNR[idx1: idx2])
            )
            lSNR.append(
                np.sum(SNR[idx1: idx2]**2)**0.5
            )
            lSNRpr.append(
                np.sum(SNR_phase_rel[idx1: idx2]**2)**0.5
            )
            if emphasis is not None:
                lemph.append(
                    np.sum(emphasis[idx1: idx2]**2)**0.5
                )

        F_Hz          = np.asarray(lF_Hz)
        data          = np.asarray(ldata)
        SNR_prev      = SNR
        SNR           = np.asarray(lSNR)
        SNR_phase_rel = np.asarray(lSNRpr)
        SNR = SNR * (np.mean(SNR_prev**2) / np.mean(SNR**2))**0.5
        if emphasis is not None:
            emphasis = np.asarray(lemph)
    ###################################### downsampling #######################

    inverse_data = grab_kwargs(aid, kw, arguments.adjustments.kw_hints, 'inverse', kwput = kwput)
    if inverse_data:
        data = 1/data

    F_boost_Hz = grab_kwargs(aid, kw, arguments.adjustments.kw_hints, 'F_boost_Hz', kwput = kwput)
    if F_boost_Hz is not None:
        if emphasis is None:
            emphasis = np.ones_like(SNR)
        for F_s, F_e in F_boost_Hz:
            select_Fb = (F_Hz >= F_s) & (F_Hz < F_e)
            emphasis[select_Fb] *= 2

    grab_kwarg_hints(aid, kw, arguments.adjustments.kw_hints, kwput = kwput)
    grab_kwarg_hints(aid, kw, arguments.ranges.kw_hints, kwput = kwput)

    ZPKrep = representations.ZPKwData(
        F_Hz          = F_Hz,
        data          = data,
        W             = SNR,
        zeros         = zeros,
        poles         = poles,
        gain          = gain,
        zeros_overlay = zeros_overlay,
        poles_overlay = poles_overlay,
        F_nyquist_Hz  = F_nyquist_Hz,
        residuals_log_im_scale = SNR_phase_rel
    )

    if aid.hint('prune_Qrank') is not None:
        aid.fitter_update(
            fitters_ZPK.ZPKrep2MRF(
                ZPKrep,
                residuals_type      = aid.hint('residuals_type'),
                coding_map          = fitters_ZPK.coding_maps.RI,
                distance_limit_auto = 2,
                delay_s             = aid.hint('delay_s'),
                delay_s_min         = aid.hint('delay_s'),
                delay_s_max         = aid.hint('delay_s'),
                h_infinity          = aid.hint('h_infinity'),
                h_infinity_deweight = aid.hint('h_infinity_deweight'),
                max_BW_Hz           = aid.hint('root_bandwidth_Hz_max'),
                F_cutoff_Hz         = aid.hint('root_F_Hz_max'),
            ),
            validate = False,
        )
        with aid.log_heading('prune Q-ranked order reduction'):
            order_reduce.order_reduce(
                aid = aid,
                Q_rank_cutoff = aid.hint('prune_Qrank'),
                reduce_c = True,
                reduce_r = True,
                optimize = False,
            )
            ZPKrep = aid.fitter.ZPKrep

    aid.fitter_update(
        fitters_ZPK.ZPKrep2MRF(
            ZPKrep,
            residuals_type      = aid.hint('residuals_type'),
            coding_map          = aid.hint('coding_map'),
            distance_limit_auto = 2,
            delay_s             = aid.hint('delay_s'),
            delay_s_min         = aid.hint('delay_s'),
            delay_s_max         = aid.hint('delay_s'),
            h_infinity          = aid.hint('h_infinity'),
            h_infinity_deweight = aid.hint('h_infinity_deweight'),
            max_BW_Hz           = aid.hint('root_bandwidth_Hz_max'),
            F_cutoff_Hz         = aid.hint('root_F_Hz_max'),
        ),
        validate = False,
    )

    if kw:
        arguments.check_remaining_arguments(kw, arguments.kw_hints)

    #TODO, put relative degree check in if the USER is specifying a filter and
    #rational fitting will not be employed

    if mode == 'dumpargs_full':
        settings                  = dict(aid.hints)
        settings['zeros']         = zeros
        settings['poles']         = poles
        settings['gain']          = gain
        settings['zeros_overlay'] = zeros_overlay
        settings['poles_overlay'] = poles_overlay
        settings['F_nyquist_Hz']  = F_nyquist_Hz
        return settings
    elif mode == 'dumpargs':
        return kwput
    elif mode == 'full':
        baseline_order = fit_full(aid, emphasis)
    elif mode == 'full2x':
        baseline_order = fit_full2x(aid, emphasis)
    elif mode == 'chebydebug':
        with aid.factorization():
            rational_fits.fit_cheby(aid)
        baseline_order = 100
    elif mode == 'chebydebug+':
        with aid.factorization():
            rational_fits.fit_cheby(aid)
        order_reduce.order_reduce(
            aid = aid,
            Q_rank_cutoff = .2,
            optimize = False,
            reduce_c = True,
            reduce_r = True,
        )
        aid.fitter.optimize(aid = aid)
        aid.fitter_update(representative = True)
        baseline_order = 100
    elif mode == 'discdebug':
        rational_fits.fit_disc(aid)
        baseline_order = 100
    elif mode == 'rational':
        baseline_order = fit_rational(aid, emphasis)
    elif mode == 'rational2x':
        baseline_order = fit_rational2x(aid, emphasis)
    elif mode == 'reduce':
        baseline_order = fit_reduce(aid, emphasis)
    elif mode == 'fit':
        baseline_order = fit_only(aid, emphasis)
    elif mode == 'copy':
        baseline_order = fit_copy(aid, emphasis)
    elif mode == 'gain':
        baseline_order = fit_copy(aid, emphasis)
        with aid.fitter.with_codings_only([aid.fitter.gain_coding]):
            aid.fitter.optimize()
        aid.fitter_update(representative = True)

    resaid = results_aid_adv.ResultsAidAdv(aid, kw = kwput)
    resaid.choose(baseline_order)
    assert(resaid.fitter is not None)
    ptbl = resaid.investigate_order_console(print_function = None)
    with aid.log_heading('investigations'):
        aid.log_info(2, ptbl)
    return resaid


def fit_full(aid, emphasis):
    SNR_adjustments.SNR_fix(aid)

    _fit_rational(aid, emphasis)

    return _reduce(aid)

def fit_full2x(aid, emphasis):
    SNR_adjustments.SNR_fix(aid)

    _fit_rational(aid, emphasis, _phase_patch = False, order_hint = 'order_first')
    _reduce(aid, with_successive = False)
    _fit_rational(aid, emphasis, _phase_patch = True)
    return _reduce(aid)


def _fit_rational(aid, emphasis, _phase_patch = True, order_hint = None):
    if emphasis is None:
        def fit_call():
            rational_fits.fit_cheby(aid, order_hint = order_hint)
            aid.log_progress(3, "Initial Order: (Z={0}, P={1}, Z-P={2})".format(
                len(aid.fitter.zeros),
                len(aid.fitter.poles),
                aid.fitter.order_relative,
            ))

            order_reduce.order_reduce(
                aid = aid,
                Q_rank_cutoff = .4,
                optimize = False,
            )

            aid.log_progress(3, "Fastdrop Order: (Z={0}, P={1}, Z-P={2})".format(
                len(aid.fitter.zeros),
                len(aid.fitter.poles),
                aid.fitter.order_relative,
            ))
            phase_patch.root_stabilize(aid)

            if _phase_patch:
                aid.log_progress(4, 'mag fitting and phase patching')
                aid.invalidate_fitters()
                phase_patch.phase_patch(aid)
                aid.fitter_update(representative = True)
    else:
        #not sure I like this emphasis mechanism
        def fit_call():
            rational_fits.fit_cheby(aid, order_hint = order_hint)
            aid.log_progress(3, "Initial Order: (Z={0}, P={1}, Z-P={2})".format(
                len(aid.fitter.zeros),
                len(aid.fitter.poles),
                aid.fitter.order_relative,
            ))

            order_reduce.order_reduce(
                aid = aid,
                Q_rank_cutoff = .4,
                optimize = False,
            )

            aid.log_progress(3, "Fastdrop Order: (Z={0}, P={1}, Z-P={2})".format(
                len(aid.fitter.zeros),
                len(aid.fitter.poles),
                aid.fitter.order_relative,
            ))

            with aid.log_heading('Applying emphasis'):
                aid.log_progress(4, 'adjusting weights')
                aid.fitter.W = aid.fitter.W * emphasis
                aid.invalidate_fitters()

                aid.fitter.optimize(aid = aid)
                aid.fitter_update(representative = False)

                algorithms.optimize_anneal(aid)
                aid.fitter_checkup()
                aid.log_progress(4, 'rational fitting (more)')
                aid.fitter_update(representative = False)

                #use the previous fits just as a suggestion
                rational_fits.fit_cheby(aid)

                order_reduce.order_reduce(
                    aid = aid,
                    Q_rank_cutoff = .4,
                    optimize = False,
                )

            if _phase_patch:
                aid.invalidate_fitters()
                phase_patch.root_stabilize(aid)
                aid.log_progress(4, 'mag fitting and phase patching')
                phase_patch.phase_patch(aid)
                aid.fitter_update(representative = True)

    with aid.log_heading('rational fitting'):
        if aid.hint('suggest'):
            #use the existing filter only as a suggestion
            fit_call()
        else:
            #apply the rational fit to the factorized form and stabilize only
            #the new cheby fit roots, not the original filter
            with aid.factorization():
                fit_call()

    aid.fitter.optimize(aid = aid)
    aid.fitter_update(representative = True)

    algorithms.optimize_anneal(aid)
    aid.fitter_checkup()
    return


def fit_rational(aid, emphasis, _phase_patch = True):
    SNR_adjustments.SNR_fix(aid)

    _fit_rational(aid, emphasis, _phase_patch = _phase_patch)

    with aid.log_heading('Q-ranked order reduction'):
        order_reduce.order_reduce(
            aid = aid,
            Q_rank_cutoff = .7,
        )
        aid.log_progress(4, 'order reduced annealing')
        algorithms.optimize_anneal(aid)
        aid.fitter_checkup()

    return aid.fitter_orders().maxzp

def fit_rational2x(aid, emphasis, _phase_patch = True):
    SNR_adjustments.SNR_fix(aid)

    _fit_rational(aid, emphasis, _phase_patch = False, order_hint = 'order_first')
    _reduce(aid, with_successive = False)
    _fit_rational(aid, emphasis, _phase_patch = _phase_patch)

    with aid.log_heading('Q-ranked order reduction'):
        order_reduce.order_reduce(
            aid = aid,
            Q_rank_cutoff = .7,
        )
        aid.log_progress(4, 'order reduced annealing')
        algorithms.optimize_anneal(aid)
        aid.fitter_checkup()

    return aid.fitter_orders().maxzp


def fit_reduce(aid, emphasis):
    SNR_adjustments.SNR_fix(aid)

    if emphasis is not None:
        aid.fitter.W = aid.fitter.W * emphasis
        aid.invalidate_fitters()

        aid.fitter.optimize(aid = aid)
        aid.fitter_update(representative = True)

        algorithms.optimize_anneal(aid)
        aid.fitter_checkup()

    aid.fitter.optimize(aid = aid)
    aid.fitter_update(representative = True)

    algorithms.optimize_anneal(aid)
    aid.fitter_checkup()

    return _reduce(aid)


def fit_only(aid, emphasis):
    SNR_adjustments.SNR_fix(aid)

    if emphasis is not None:
        aid.fitter.W = aid.fitter.W * emphasis
        aid.invalidate_fitters()
        with aid.fitter.with_codings_only([aid.fitter.gain_coding]):
            aid.fitter.optimize()
        aid.fitter_update(representative = True)

        #algorithms.optimize_anneal(aid)
        aid.fitter.optimize(aid = aid)
        aid.fitter_update(representative = True)

        algorithms.optimize_anneal(aid)
        aid.fitter_checkup()
    else:
        with aid.fitter.with_codings_only([aid.fitter.gain_coding]):
            aid.fitter.optimize()
        aid.fitter_update(representative = True)

    aid.fitter.optimize(aid = aid)
    aid.fitter_update(representative = True)

    algorithms.optimize_anneal(aid)
    aid.fitter_checkup()

    if aid.hint('delay_s_max') is not None:
        aid.log_progress(4, 'activating delay fitting')
        #print('delay_coding: ', aid.fitter.delay_s_max)
        #print('delay_coding: ', aid.fitter.codings_revision)
        aid.fitter.delay_s_max = aid.hint('delay_s_max')
        aid.fitter.delay_s_min = aid.hint('delay_s_min')
        algorithms.optimize_anneal(aid)
        #print('delay_coding: ', aid.fitter.codings_revision)
        #print('delay_coding: ', aid.fitter.delay_s_max)
        #print('delay_coding: ', aid.fitter.delay_coding.N_parameters)
        aid.fitter_checkup()

    algorithms.optimize_anneal(aid)
    aid.fitter_checkup()

    if aid.hint('delay_s_max') is not None:
        aid.log_alert(2, "Baseline fit delay: ", aid.fitter.delay_s)
    baseline_order = aid.fitter_orders().maxzp
    aid.log_alert(2, (
        "Baseline fit residuals: {:.2e}, at order {}"
    ).format(aid.fitter.residuals_average, baseline_order))
    return baseline_order


def fit_copy(aid, emphasis):
    SNR_adjustments.SNR_fix(aid)
    aid.fitter_update(representative = True)
    baseline_order = aid.fitter_orders().maxzp
    aid.log_alert(2, (
        "Baseline fit residuals: {:.2e}, at order {}"
    ).format(aid.fitter.residuals_average, baseline_order))
    return baseline_order


def _reduce(aid, with_successive = True):
    with aid.log_heading('Q-ranked order reduction'):
        rzp_list = order_reduce.order_reduce(
            aid = aid,
            Q_rank_cutoff = 2,
        )
        aid.log_progress(4, 'order reduced annealing')
        algorithms.optimize_anneal(aid)
        improved = aid.fitter_checkup()

    with aid.log_heading('Q-ranked selective order restoration'):
        if improved:
            rzp_list.sort()
            order_reduce.order_restore(aid, rzp_list)

    algorithms.optimize_anneal(aid)
    aid.fitter_checkup()

    with aid.log_heading('selective order reduction'):
        order_reduce.order_reduce_selective(
            aid = aid,
            num_total_max = 8,
            num_type_max  = 2,
        )

    if aid.hint('delay_s_max') is not None:
        aid.log_progress(4, 'activating delay fitting')
        aid.fitter.delay_s_max = aid.hint('delay_s_max')
        aid.fitter.delay_s_min = aid.hint('delay_s_min')
        algorithms.optimize_anneal(aid)
        aid.fitter_checkup()

    algorithms.optimize_anneal(aid)
    aid.fitter_checkup()

    if aid.hint('delay_s_max') is not None:
        aid.log_alert(2, "Baseline fit delay: ", aid.fitter.delay_s)
    baseline_order = aid.fitter_orders().maxzp
    aid.log_alert(2, (
        "Baseline fit residuals: {:.2e}, at order {}"
    ).format(aid.fitter.residuals_average, baseline_order))

    if aid.hint('baseline_only'):
        return baseline_order
    #TODO, now need to make a choice based on residuals
    if with_successive:
        with aid.log_heading('successive order reduction'):
            order_reduce.order_reduce_successive(
                aid = aid,
                num_total_max = 6,
                num_type_max  = 2,
            )
    return baseline_order


data2filter.__doc__ = autogen_docstr.__doc__
