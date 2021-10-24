# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest
import numpy as np

import os.path as path

from wavestate.iirrational import v2
from wavestate.iirrational.v2 import testing

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
def test_simple1(request, browser, plotsections, plot_verbosity):
    data_name = "PUM_long"
    data = IIRrational_data_dev(
        data_name,
        instance_num=1,
        set_num=1,
    )
    out = v2.data2filter(
        xfer=data.data,
        F_Hz=data.F_Hz,
        SNR=data.SNR,
        order_initial=10,
        delay_s=0,
        F_nyquist_Hz=None,
        poles=(-100, -50),
        zeros=(-10, -50),
        mode="fit",
        delay_s_max=0.01,
        h_infinity=0.99,
        h_infinity_deweight=0.1,
        # hints = [
        #    testing.validate_plot_log(__file__, request),
        # ],
    )

    print(out.as_scipy_signal_ZPKsw())
    print(out.as_foton_ZPKsf())
    print(out.as_foton_str_ZPKsf())
    print(out.as_foton_ZPKsn())
    print(out.as_foton_str_ZPKsn())
    return


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_simple1_vec(request, browser, plotsections, plot_verbosity):
    data_name = "PUM_long"
    data = IIRrational_data_dev(
        data_name,
        instance_num=1,
        set_num=1,
    )
    out = v2.data2filter(
        xfer=data.data,
        F_Hz=data.F_Hz,
        SNR=data.SNR,
        order_initial=10,
        delay_s=0,
        F_nyquist_Hz=None,
        poles=(-100, -50),
        zeros=(-10, -50),
        mode="fit",
        delay_s_max=0.01,
        h_infinity=abs(1 / (1 + 1j * data.F_Hz / 1)),
        # hints = [
        #    testing.validate_plot_log(__file__, request),
        # ],
    )

    print(out.as_scipy_signal_ZPKsw())
    print(out.as_foton_ZPKsf())
    print(out.as_foton_str_ZPKsf())
    print(out.as_foton_ZPKsn())
    print(out.as_foton_str_ZPKsn())
    return


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_simple1_callable(request, browser, plotsections, plot_verbosity):
    data_name = "PUM_long"
    data = IIRrational_data_dev(
        data_name,
        instance_num=1,
        set_num=1,
    )
    out = v2.data2filter(
        xfer=data.data,
        F_Hz=data.F_Hz,
        SNR=data.SNR,
        order_initial=10,
        delay_s=0,
        F_nyquist_Hz=None,
        poles=(-100, -50),
        zeros=(-10, -50),
        mode="fit",
        delay_s_max=0.01,
        h_infinity=lambda L: abs(1 / (1 + 1j * np.linspace(0, 100, L))),
        # hints = [
        #    testing.validate_plot_log(__file__, request),
        # ],
    )

    print(out.as_scipy_signal_ZPKsw())
    print(out.as_foton_ZPKsf())
    print(out.as_foton_str_ZPKsf())
    print(out.as_foton_ZPKsn())
    print(out.as_foton_str_ZPKsn())
    return
