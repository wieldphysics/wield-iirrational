# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

# from wield.iirrational import fitters
from wield.iirrational import fitters_rational
from wield.iirrational import plots
from wield.iirrational.utilities import relpy
import numpy as np
from os import path

try:
    from IIRrational_test_data import IIRrational_data_dev
    from IIRrational_test_data import matlab as matlab_test_data
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False

localpath = path.split(__file__)[0]
skipdeco = pytest.mark.skipif(
    module_import_skip, reason="cannot import IIRrational_test_data"
)


@skipdeco
def test_imcheby_convergence():
    data = IIRrational_data_dev("quad_M0_P2V_10Hz")
    SNR = np.minimum(data.SNR, 100)
    fitter = fitters_rational.ChebychevFilter(
        F_Hz=data.F_Hz,
        data=data.data,
        W=SNR,
        npoles=100,
        nzeros=100,
    )
    print("data", fitter.data)
    print("W", fitter.W)
    print("h_a", fitter.h_a.xfer, fitter.h_a.lnG)
    print("h_b", fitter.h_b.xfer, fitter.h_b.lnG)
    fitter.fit_zeros()
    fitter.fit_poles()
    print("h_a", fitter.h_a.xfer[0], fitter.h_a.lnG)
    print("h_b", fitter.h_b.xfer[0], fitter.h_b.lnG)
    fitter.fit_zeros()
    print("h_a", fitter.h_a.xfer[0], fitter.h_a.lnG)
    print("h_b", fitter.h_b.xfer[0], fitter.h_b.lnG)
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff=0.3)
    fitter.stabilize = True
    # fitter.npoles = len(fitter.poles)
    # fitter.nzeros = len(fitter.zeros)
    ##fitter.stabilize = True
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff=0.7)

    fitter2 = fitters_rational.ChebychevFilter(
        F_Hz=fitter.F_Hz,
        data=fitter.data,
        W=fitter.W,
    )
    # fitter2.a_vec = fitter.a_vec
    # fitter2.b_vec = fitter.b_vec
    # fitter2.gain = fitter.gain
    # fitter2.poles = fitter.poles
    # fitter2.zeros = fitter.zeros
    # fitter2.gain = fitter.gain

    # TODO
    #!Make this a TEST so it doesn't plot all of the time

    axB = plots.plot_fit(
        fitter,
        plot_zp=True,
        # plot_past_data = False,
        plot_past_data=True,
        xscale="log",
    )

    axB.save(relpy(__file__, "plots/imcheby-convergence"))
    return


@skipdeco
def test_imcheby_convergence_2():
    data = IIRrational_data_dev("quad_M0_P2V_10Hz")
    SNR = np.minimum(data.SNR, 100)
    fitter = fitters_rational.ChebychevFilter(
        F_Hz=data.F_Hz,
        data=data.data,
        W=SNR,
        npoles=20,
        nzeros=18,
    )
    print("data", fitter.data)
    print("W", fitter.W)
    print("h_a", fitter.h_a.xfer, fitter.h_a.lnG)
    print("h_b", fitter.h_b.xfer, fitter.h_b.lnG)
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    ##print("HOLD")
    fitter.npoles = 62
    fitter.nzeros = 60
    ####print('h_a', fitter.h_a.xfer[0], fitter.h_a.lnG)
    ####print('h_b', fitter.h_b.xfer[0], fitter.h_b.lnG)
    fitter.fit_zeros()
    ####print('h_a', fitter.h_a.xfer[0], fitter.h_a.lnG)
    ####print('h_b', fitter.h_b.xfer[0], fitter.h_b.lnG)
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.matched_pairs_clear(Q_rank_cutoff=0.3)
    fitter.stabilize = True
    ####fitter.npoles = len(fitter.poles)
    ####fitter.nzeros = len(fitter.zeros)
    #####fitter.stabilize = True

    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff=0.7)
    # fitter.fit_poles()
    # fitter.fit_zeros()

    # TODO
    #!Make this a TEST so it doesn't plot all of the time

    axB = plots.plot_fit(
        fitter,
        plot_zp=True,
        plot_past_data=False,
        # plot_past_data = True,
        xscale="log",
    )
    axB.save(relpy(__file__, "plots/imcheby-convergence_2"))

    axB = plots.plot_fit(
        fitter,
        plot_zp=True,
        # plot_past_data = False,
        plot_past_data=True,
        xscale="log",
    )

    axB.save(relpy(__file__, "plots/imcheby-convergence_3"))
    return


@skipdeco
def test_imcheby_convergence_alpha():
    fname = matlab_test_data.matfiles["alphas_fdomain.txt"]
    alphas = np.loadtxt(fname)
    # assemble back as complex numbers (each column of a is a TF)
    fr = alphas[:, 0]
    a = np.zeros((alphas.shape[0], (alphas.shape[1] - 1) // 2), dtype=complex)
    for i in range(a.shape[1]):
        a[:, i] = alphas[:, i + 1] + 1j * alphas[:, i + 2]

    idx = np.where(np.logical_and(fr > 30, fr < 300))[0]

    order = 20
    fitter = fitters_rational.ChebychevFilter(
        data=a[idx, 0],
        F_Hz=fr[idx],
        W=10,
        zeros_overlay=(-1,),
        poles_overlay=(-30,),
        npoles=order // 2,
        nzeros=order // 2,
    )
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.nzeros = order
    fitter.npoles = order
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff=0.5)
    fitter.stabilize = True
    # fitter.npoles = len(fitter.poles)
    # fitter.nzeros = len(fitter.zeros)
    ##fitter.stabilize = True
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_SVD()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff=0.5)
    # fitter.stabilize = False
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.fit_poles()
    fitter.fit_zeros()
    fitter.matched_pairs_clear(Q_rank_cutoff=0.1)

    print("Z: ", fitter.ZPKrep.zeros)
    print("P: ", fitter.ZPKrep.poles)

    # from wield.iirrational import fitters_ZPK
    # fitter = fitters_ZPK.ZPKrep2MRF(
    #    ZPKrep = fitter.ZPKrep,
    #    residuals_type      = 'log',
    #    coding_map          = fitters_ZPK.codings_s.coding_maps.SOSnl,
    #    #coding_map          = fitters_ZPK.codings_s.coding_maps.nlFBW,
    #    distance_limit_auto = 0,
    #    max_BW_Hz           = np.max(fitter.F_Hz),
    # )
    # TODO
    #!Make this a TEST so it doesn't plot all of the time

    axB = plots.plot_fit(
        fitter,
        plot_zp=True,
        # plot_past_data = False,
        plot_past_data=True,
        xscale="log",
    )

    axB.save(relpy(__file__, "plots/imcheby-convergence"))
    return
