# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os.path as path

from IIRrational.testing import IIRrational_data
from IIRrational import v2
from IIRrational.v2 import testing

from IIRrational.testing.utilities import (
    sign_validate_and_plot_hint,
)

def test_simple1(request, browser, plotsections, plot_verbosity):
    data_name = 'simple1'
    data = IIRrational_data(
        data_name,
        instance_num = 1,
        set_num = 1,
    )
    out = v2.data2filter(
        xfer          = data.data,
        F_Hz          = data.F_Hz,
        SNR           = data.SNR,
        order_initial = 10,
        delay_s       = 0,
        F_nyquist_Hz  = None,
        poles = (-1, -1, -1, -1, -1),
        zeros = (-.1, -.1, -.1),
        mode = 'fit',
        #hints = [
        #    testing.validate_plot_log(__file__, request),
        #],
    )

    print(out.as_scipy_signal_ZPKsw())
    print(out.as_foton_ZPKsf())
    print(out.as_foton_str_ZPKsf())
    print(out.as_foton_ZPKsn())
    print(out.as_foton_str_ZPKsn())
    return
