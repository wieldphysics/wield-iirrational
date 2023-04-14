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

try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_disc_convergence():
    data = IIRrational_data_dev("quad_M0_P2V_10Hz")
    SNR = np.minimum(data.SNR, 100)
    F_nyquist_Hz = np.max(data.F_Hz)
    fitter = fitters_rational.RationalDiscFilter(
        F_Hz=data.F_Hz,
        data=data.data,
        W=SNR,
        F_nyquist_Hz=F_nyquist_Hz,
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

    # TODO
    #!Make this a TEST so it doesn't plot all of the time

    axB = plots.plot_fit(
        fitter,
        plot_zp=True,
        # plot_past_data = False,
        plot_past_data=True,
        xscale="log",
    )

    axB.save(relpy(__file__, "plots/disc-convergence"))
    return
