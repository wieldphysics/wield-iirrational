# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
from IIRrational_test_data import IIRrational_data_dev
from IIRrational.utilities.mpl import mplfigB
from IIRrational import testing

from IIRrational.fitters_ZPK import ZPKrep2MRF
from IIRrational.fitters_ZPK import codings_s
from IIRrational.fitters_ZPK.MRF import check_jac

from utilities import check_residuals


def prep_fitter():
    data = IIRrational_data_dev("quad_M0_P2V_10Hz")

    from IIRrational import fitters_rational

    SNR = np.minimum(data.SNR, 100)
    fitter = fitters_rational.ChebychevFilter(
        F_Hz=data.F_Hz,
        data=data.data,
        W=SNR,
        npoles=100,
        nzeros=100,
    )
    # fitter.poles = (.90,) * 0
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_SVD()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    # fitter.matched_pairs_clear()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.stabilize = True
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.matched_pairs_clear()
    # axB = plot_fit(fitter, with_error=False, plot_zp = True)
    # axB = mplfigB()
    # axB.ax0.semilogy(fitter.F_Hz, abs(fitter.xfer_fit))
    return fitter


def test_NL_fitter1():
    fitter = prep_fitter()
    rep_s = ZPKrep2MRF(
        fitter.ZPKrep,
        coding_map=codings_s.coding_maps.nlFBW_safe,
    )

    print(rep_s.residuals_average)
    # axB = mplfigB()

    check_residuals(rep_s)
    opt = rep_s.optimize()
    print(rep_s.residuals_average)
    assert rep_s.residuals_average < 200

    # axB = plot_fitter_flag(rep_s, with_error=False, plot_zp = True)
