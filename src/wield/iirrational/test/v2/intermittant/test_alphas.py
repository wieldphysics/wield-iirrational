# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest
import numpy as np
from os import path

import IIRrational
from wavestate.iirrational import v2
from wavestate.iirrational.testing.plots import plot_on_assert

from wavestate.iirrational.testing.utilities import (
    sign_validate_and_plot_hint,
    stability_validate_and_plot_hint,
)

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
def test_alphas1(request, browser, plotsections, plot_verbosity):
    fname = matlab_test_data.matfiles["alphas_fdomain.txt"]
    alphas = np.loadtxt(fname)
    # assemble back as complex numbers (each column of a is a TF)
    fr = alphas[:, 0]
    a = np.zeros((alphas.shape[0], (alphas.shape[1] - 1) // 2), dtype=complex)
    for i in range(a.shape[1]):
        a[:, i] = alphas[:, i + 1] + 1j * alphas[:, i + 2]

    idx = np.where(np.logical_and(fr > 30, fr < 300))[0]

    fit = wavestate.iirrational.v2.data2filter(
        data=a[idx, 0],
        F_Hz=fr[idx],
        SNR=100,
        # mode = 'rational',
        mode="chebydebug",
        F_nyquist_Hz=None,
        log_level=10,
        baseline_only=True,
        # poles_overlay = (-1,),
        # zeros_overlay = (-30,),
        # never_unstable_poles = True,
        # never_unstable_zeros = True,
        log_level_debug=10,
    )
    with plot_on_assert(__file__, request, fit.fitter, plot_anyway=True):
        pass
    return


@skipdeco
def test_alphas1_v1good(request, browser, plotsections, plot_verbosity):
    fname = matlab_test_data.matfiles["alphas_fdomain.txt"]
    alphas = np.loadtxt(fname)
    # assemble back as complex numbers (each column of a is a TF)
    fr = alphas[:, 0]
    a = np.zeros((alphas.shape[0], (alphas.shape[1] - 1) // 2), dtype=complex)
    for i in range(a.shape[1]):
        a[:, i] = alphas[:, i + 1] + 1j * alphas[:, i + 2]

    idx = np.where(np.logical_and(fr > 30, fr < 300))[0]

    fit = wavestate.iirrational.v1.data2filter(
        data=a[idx, 0],
        F_Hz=fr[idx],
        SNR=100,
        order_initial=50,
        F_nyquist_Hz=None,
    )
    with plot_on_assert(__file__, request, fit.fitter, plot_anyway=True):
        pass
    return
