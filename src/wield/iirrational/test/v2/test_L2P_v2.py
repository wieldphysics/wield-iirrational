# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from wavestate.iirrational import v2
from wavestate.iirrational.testing.plots import plot_on_assert

from wavestate.iirrational.testing.utilities import (
    sign_validate_and_plot_hint,
)

try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_L2P_MULT_S(request, browser, plotsections, plot_verbosity):
    data_name = "L2P_MULT_WLF"
    dataB = IIRrational_data_dev(
        data_name,
        instance_num=1,
        set_num=2,
    )
    out = v2.data2filter(
        F_Hz=dataB.F_Hz,
        data=dataB.data,
        SNR=dataB.SNR,
        F_nyquist_Hz=None,
        order_initial=15,
        # order_initial = 80,
        baseline_only=True,
        delay_s=0,
        hints=[sign_validate_and_plot_hint(__file__, request)],
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway=True):
        pass
    return


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_L2P_MULT_S_WOLAY(request, browser, plotsections, plot_verbosity):
    data_name = "L2P_MULT_WLF"
    dataB = IIRrational_data_dev(
        data_name,
        instance_num=1,
        set_num=2,
    )
    out = v2.data2filter(
        F_Hz=dataB.F_Hz,
        data=dataB.data,
        SNR=dataB.SNR,
        F_nyquist_Hz=None,
        order_initial=30,
        delay_s=0,
        ZPK_overlay=((0,), (), 1),
        baseline_only=True,
        hints=[sign_validate_and_plot_hint(__file__, request)],
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway=True):
        print(out.fitter)
        print(out.fitter.F_nyquist_Hz)
        print(out.fitter.zeros)
        print(out.fitter.zeros_overlay)
        # print(out.fitter.ZPKrep)
        # assert(False)
        pass
    print(out.as_foton_str_ZPKsf(annotate_pairs=False))
    return
