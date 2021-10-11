# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
import declarative

from .. import fitters_rational
from .. import annotate
from .. import plots
from ..utilities import args, ensure_aid

from . import nyquist_move


def rational_disc_fit_mag(
    argB             = args.UNSPEC,
    F_Hz             = args.UNSPEC,
    data             = args.UNSPEC,
    nyquist_final_Hz = args.UNSPEC,
    SNR              = args.UNSPEC,
    SNR_cutoff       = args.UNSPEC,
    order            = args.UNSPEC,
    order_max        = args.UNSPEC,
    doc_db           = None,
    aid              = None,
):
    aid              = ensure_aid(aid)

    F_Hz             = args.argscan(locals(), argB, args.REQ, arg = 'F_Hz')
    data             = args.argscan(locals(), argB, args.REQ, arg = 'data')
    nyquist_final_Hz = args.argscan(locals(), argB, None,     arg = 'nyquist_final_Hz')
    SNR              = args.argscan(locals(), argB, args.REQ, arg = 'SNR')
    SNR_cutoff       = args.argscan(locals(), argB, 100,      arg = 'SNR_cutoff')
    order            = args.argscan(locals(), argB, None,     arg = 'order')
    order_max        = args.argscan(locals(), argB, 100 ,     arg = 'order_max')

    if order is not None:
        order     = int(order)

    if order_max is not None:
        order_max = int(order_max)

    N_first = 20
    N_final = int(min(len(F_Hz) // 10, order_max))

    if order is not None:
        return SVD_method_mag(
            F_Hz             = F_Hz,
            data             = data,
            nyquist_final_Hz = nyquist_final_Hz,
            SNR              = SNR,
            SNR_cutoff       = SNR_cutoff,
            order            = order,
            doc_db           = doc_db,
            aid              = aid,
        )

    if doc_db is None:
        doc_db_last = None
        doc_db_current = None

    doc_db_true = doc_db
    #otherwise, scan through
    N_current = N_first

    if doc_db is not None:
        doc_db_last = declarative.DeepBunch()
        doc_db_last.update_recursive(doc_db)
    fitter_last = SVD_method_mag(
        F_Hz             = F_Hz,
        data             = data,
        nyquist_final_Hz = None,
        SNR              = SNR,
        SNR_cutoff       = SNR_cutoff,
        order            = N_current ,
        doc_db           = doc_db_last,
        aid              = aid,
    )
    restot_last = fitter_last.residuals_average

    def fitter_ord(fitter):
        return max(len(fitter.zeros), len(fitter.poles))

    while True:
        if N_current == N_final:
            fitter_use = fitter_last
            break

        N_current = N_current * 2
        if N_current > N_final:
            N_current = N_final

        if doc_db is not None:
            doc_db_current = declarative.DeepBunch()
            doc_db_current.update_recursive(doc_db)
        fitter = SVD_method_mag(
            F_Hz             = F_Hz      ,
            data             = data      ,
            nyquist_final_Hz = None      ,
            SNR              = SNR       ,
            SNR_cutoff       = SNR_cutoff,
            order            = N_current ,
            doc_db           = doc_db_current,
            aid              = aid,
        )
        restot = fitter.residuals_average
        fitter_red = fitter.copy()
        fitter_red.matched_pairs_clear(Q_rank_cutoff = .2)
        restot_red = fitter_red.residuals_average

        if restot_last < restot:
            fitter_use = fitter_last
            if doc_db is not None:
                doc_db.update_recursive(doc_db_last)
            aid.log("Using last (direct)!", fitter_ord(fitter_last))
            break

        if restot_last < restot_red:
            fitter_use = fitter_last
            if doc_db is not None:
                doc_db.update_recursive(doc_db_last)
            aid.log("Using last (reduced)!", fitter_ord(fitter_last))
            break

        if restot_last < 1.10 * restot:
            ord_red = fitter_ord(fitter_red)
            ord_last = fitter_ord(fitter_last)
            if ord_red < ord_last:
                fitter_use = fitter_red
                if doc_db is not None:
                    doc_db.update_recursive(doc_db_current)
            else:
                fitter_use = fitter_last
                if doc_db is not None:
                    doc_db.update_recursive(doc_db_last)
            aid.log("Using current")
            break
        #else continue the loop
        fitter_last = fitter
        restot_last = restot
        doc_db_last = doc_db_current

    fitter = fitter_use
    if nyquist_final_Hz is not None:
        nyquist_move.nyquist_move(
            fitter,
            nyquist_new = nyquist_final_Hz,
            split_neg_real = True,
            clear_neg_real = True,
            aid              = aid,
        )
    return fitter


def SVD_method_mag(
    argB             = args.UNSPEC,
    F_Hz             = args.UNSPEC,
    data             = args.UNSPEC,
    nyquist_final_Hz = args.UNSPEC,
    SNR              = args.UNSPEC,
    SNR_cutoff       = args.UNSPEC,
    order            = args.UNSPEC,
    doc_db           = None,
    aid              = None,
):
    aid = ensure_aid(aid)

    F_Hz             = args.argscan(locals(), argB, args.REQ, arg = 'F_Hz')
    data             = args.argscan(locals(), argB, args.REQ, arg = 'data')
    nyquist_final_Hz = args.argscan(locals(), argB, None,     arg = 'nyquist_final_Hz')
    SNR              = args.argscan(locals(), argB, args.REQ, arg = 'SNR')
    SNR_cutoff       = args.argscan(locals(), argB, 100,      arg = 'SNR_cutoff')
    order            = args.argscan(locals(), argB, 50,       arg = 'order')

    F_Hz_sort    = np.sort(F_Hz)
    F_nyquist_Hz = (F_Hz_sort[-1] + (F_Hz_sort[-1] - F_Hz_sort[-2]))
    F_nyquist_Hz = 1 * F_Hz_sort[-1]

    if SNR_cutoff is not None:
        if SNR is None:
            SNR = 1
        else:
            SNR = np.minimum(SNR, SNR_cutoff)

    kwargs = dict(
        #delay_s = -.5 / F_nyquist_Hz
    )

    fitter_x = fitters_rational.RationalDiscFilterMag(
        F_Hz         = F_Hz,
        data         = data,
        W            = SNR,
        npoles       = order,
        nzeros       = order,
        F_nyquist_Hz = F_nyquist_Hz,
        **kwargs
    )
    fitter_x.fit_zeros()
    fitter_x.fit_poles()
    fitter_x.fit_zeros()
    fitter_x.fit_poles()
    fitter_x.fit_zeros()
    annotate.annotate(
        doc_db,
        name      = 'initial_direct',
        fitter    = fitter_x,
        plotter   = plots.plot_fitter_flag,
        method    = 'fit_poles, fit_zeros',
        verbosity = 9,
        about = (
            """
            initial guess without SVD technique
            """
        ),
    )

    fitter_p = fitters_rational.RationalDiscFilterMag(
        #parent       = fitter_x,
        F_Hz         = F_Hz,
        data         = data,
        W            = SNR,
        npoles       = order,
        nzeros       = order,
        F_nyquist_Hz = F_nyquist_Hz,
        **kwargs
    )
    fitter_p.fit_poles_mod_zeros()
    fitter_p.fit_zeros()
    fitter_p.fit_poles()
    fitter_p.fit_zeros()
    fitter_p.fit_poles()
    fitter_p.fit_zeros()
    fitter_p.fit_poles()
    fitter_p.fit_zeros()
    annotate.annotate(
        doc_db,
        name      = 'initial_poles',
        fitter    = fitter_p,
        plotter   = plots.plot_fitter_flag,
        method    = 'fit_poles_mod_zeros',
        verbosity = 9,
        about = (
            """
            Performs the SVD for a rough initial guess
            """
        ),
    )

    #
    rep_z = fitters_rational.RationalDiscFilterMag(
        #parent       = fitter_x,
        F_Hz         = F_Hz,
        data         = data,
        W            = SNR,
        npoles       = order,
        nzeros       = order,
        F_nyquist_Hz = F_nyquist_Hz,
        **kwargs
    )
    #put on mindelay for zeros to help ensure stable poles in froissart doublets
    rep_z.fit_zeros_mod_poles()
    rep_z.fit_poles()
    rep_z.fit_zeros()
    rep_z.fit_poles()
    rep_z.fit_zeros()
    rep_z.fit_poles()
    rep_z.fit_zeros()
    annotate.annotate(
        doc_db,
        name      = 'initial_zeros',
        fitter    = rep_z,
        plotter   = plots.plot_fitter_flag,
        method    = 'fit_poles_mod_zeros',
        verbosity = 9,
        about = (
            """
            Performs the SVD for a rough initial guess
            """
        ),
    )

    aid.log("(direct = {0:.3e}, Psvd= {0:.3e}, Zsvd= {0:.3e})".format(fitter_x.residuals_average, fitter_p.residuals_average, rep_z.residuals_average,))
    fitter_list = [
        fitter_x,
        #fitter_p,
        #rep_z,
    ]
    fitter_list.sort(key = lambda f : f.residuals_average)
    fitter = fitter_list[0]
    if fitter is fitter_p:
        annotate.annotate(
            doc_db,
            name      = 'choose poles',
            method    = 'if',
            verbosity = 3,
            about = (
                """
                Chose the Poles SVD fitter as it had the smaller residual of {0:.2e} vs. {1:.2e} for the zeros
                """.format(fitter_p.residuals_average, rep_z.residuals_average)
            ),
        )
        #fitter.fit_poles_mod_zeros()
        #fitter.fit_zeros()
    elif fitter is rep_z:
        annotate.annotate(
            doc_db,
            name      = 'choose zeros',
            method    = 'if',
            verbosity = 3,
            about = (
                """
                Chose the zeros SVD fitter as it had the smaller residual of {0:.2e} vs. {1:.2e} for the poles
                """.format(rep_z.residuals_average, fitter_p.residuals_average)
            ),
        )
        #fitter.fit_zeros_mod_poles()
        #fitter.fit_poles()
    elif fitter is fitter_x:
        annotate.annotate(
            doc_db,
            name      = 'choose direct',
            method    = 'if',
            verbosity = 3,
            about = (
                """
                Chose the direct (non-SVD) fitter as it had the smallest residual of {0:.2e}
                """.format(fitter_x.residuals_average)
            ),
        )
    #return fitter
    #fitter.match_pair_iter(Q_rank_cutoff = .4, zeros_first = False)
    fitter.fit_poles()
    fitter.fit_zeros()
    #fitter.match_pair_iter(Q_rank_cutoff = .4)
    fitter.fit_poles()
    fitter.fit_zeros()

    r_last = 1e10
    for idx in range(3):
        #fitter.match_pair_iter(Q_rank_cutoff = .4)
        fitter.fit_poles()
        r_p = fitter.residuals_average
        fitter.fit_zeros()
        r_z = fitter.residuals_average
        r = min(r_p, r_z)
        #aid.log(r / r_last)
        r_last = r

    annotate.annotate(
        doc_db,
        name     = 'seq_iter_3',
        fitter   = fitter,
        plotter  = plots.plot_fitter_flag,
        method   = 'RationalDiscFilterMag.fit_poles, RationalDiscFilterMag.fit_zeros',
        about    = (
            """
            First iterations, enforcing stabilized poles residual of {0:.2e}
            """.format(fitter.residuals_average)
        ),
        verbosity       = 9,
    )

    fitter.fit_zeros()
    r_z = fitter.residuals_average
    fitter.fit_poles()
    r_p = fitter.residuals_average
    if r_z < r_p:
        #FIT POLES SHOULD BE THE LAST FIT! The residuals tend to be lowest under the poles
        fitter.fit_zeros()
        pass

    annotate.annotate(
        doc_db,
        name     = 'seq_iter_4',
        fitter   = fitter,
        plotter  = plots.plot_fitter_flag,
        method   = 'RationalDiscFilterMag.fit_poles, RationalDiscFilterMag.fit_zeros',
        about    = (
            """
            First iterations, enforcing stabilized poles residual of {0:.2e}
            """.format(fitter.residuals_average)
        ),
        verbosity       = 4,
    )
    aid.log("LINEAR Final Residuals: ", fitter.residuals_average)

    if nyquist_final_Hz is not None:
        fitter = fitters_rational.RationalDiscFilter(
            F_Hz         = F_Hz,
            data         = data,
            W            = SNR,
            npoles       = order,
            nzeros       = order,
            F_nyquist_Hz = F_nyquist_Hz,
            ZPK          = fitter.ZPK,
            **kwargs
        )
        ZPK = nyquist_move.nyquist_move(
            fitter,
            nyquist_new = nyquist_final_Hz,
            split_neg_real = True,
            clear_neg_real = True,
            aid              = aid,
        )

    #annotate.annotate(
    #    doc_db,
    #    name     = 'nyquist_move',
    #    fitter   = fitter,
    #    method   = 'RationalDiscFilterMag.nyquist_move',
    #    plotter  = plots.plot_fitter_flag,
    #    about    = (
    #        """
    #        Moves to intended nyquist, clipping negative real poles/zeros likely caused by phase defects of low nyquist
    #        """.format(fitter.residuals_average)
    #    ),
    #    verbosity       = 3,
    #)

    #NOTE: shouldn't currently calculate residuals_average using the post-nyquist fitter as it is not numerically stable
    #TODO: switch to MRF for this residual
    annotate.annotate(
        doc_db,
        name     = "Final",
        method   = "SVD_method",
        fitter   = fitter,
        plotter  = plots.plot_fitter_flag,
        about    = (
            """
            Uses SVD to create initial guess of fit for data, followed by several iterative fits
            (see reference ???).

               * It is a linear method, finding global optimum (nonlocal). This makes it get stuck if systematics are bad. To prevent this,
            it requires gratuitous overfitting to reliably get good fits.
               * It requires a nyquist frequency that is very low, near the last data point. This can cause artifacts due to phasing discontinuity near the nyquist.
               * The provided nyquist frequency is shifted up at the end, removing the real poles/zeros that are typically due to phasing discontinuity

            """.format()
        ),
        verbosity       = 2,
    )
    return fitter


