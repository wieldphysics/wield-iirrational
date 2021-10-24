# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from wavestate.iirrational import v1
from wavestate.iirrational.testing.utilities import (
    sign_validate_hint,
    sign_validate_and_plot_hint,
)
from wavestate.iirrational import fitters_rational, fitters_ZPK
from wavestate.iirrational.testing.plots import plot_on_assert

try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_L2P_MULT(request, browser, plotsections, plot_verbosity):
    data_name = "L2P_MULT"
    out = v1.data2filter(
        IIRrational_data_dev(
            data_name,
            instance_num=1,
            set_num=2,
        ),
        hints=[sign_validate_and_plot_hint(__file__, request)],
        # order_initial = 80,
    )

    with plot_on_assert(__file__, request, out.fitter, plot_anyway=True):
        pass
    return


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_L2P_MULT_S(request, browser, plotsections, plot_verbosity):
    data_name = "L2P_MULT_WLF"
    out = v1.data2filter(
        IIRrational_data_dev(
            data_name,
            instance_num=1,
            set_num=2,
        ),
        F_nyquist_Hz=None,
        # order_initial = 80,
        hints=[sign_validate_and_plot_hint(__file__, request)],
        delay_s=0,
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway=True):
        pass
    return


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_L2P_MULT_S_WOLAY(request, browser, plotsections, plot_verbosity):
    data_name = "L2P_MULT_WLF"
    out = v1.data2filter(
        IIRrational_data_dev(
            data_name,
            instance_num=1,
            set_num=2,
        ),
        F_nyquist_Hz=None,
        order_initial=80,
        delay_s=0,
        ZPK_overlay=((0,), (), 1),
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
    return


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_L2P_MULT_WOLAY(request, browser, plotsections, plot_verbosity):
    data_name = "L2P_MULT_WLF"
    out = v1.data2filter(
        IIRrational_data_dev(
            data_name,
            instance_num=1,
            set_num=2,
        ),
        delay_s=0,
        ZPK_overlay=((0.99999,), (), 1),
        hints=[sign_validate_and_plot_hint(__file__, request)],
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway=True):
        print(out.fitter)
        print(out.fitter.F_nyquist_Hz)
        print(out.fitter.zeros)
        print(out.fitter.zeros_overlay)
        print(out.fitter.ZPKrep)
        # assert(False)
        pass
    return


if __name__ == "__main__":
    test_simple1(False, None, 10)
