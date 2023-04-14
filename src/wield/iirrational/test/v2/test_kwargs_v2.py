# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os.path as path

from wield.iirrational.testing import IIRrational_data
from wield.iirrational import v2
from wield.iirrational.v2.__main__ import main as v2_main
from wield.iirrational.v2 import testing

from wield.iirrational.testing.utilities import (
    sign_validate_and_plot_hint,
)


def test_argtell(request, browser, plotsections, plot_verbosity):
    dset = IIRrational_data(
        "simple1",
        instance_num=1,
        set_num=1,
    )

    res = v2.data2filter(
        F_Hz=dset.F_Hz,
        SNR=dset.SNR,
        data=dset.data,
        order_initial=10,
        delay_s=0,
        F_nyquist_Hz=None,
        mode="fit",
        delay_f=5,
        hints=[
            testing.validate_plot_log(__file__, request),
        ],
    )
    return


def test_prog_help(request, browser, plotsections, plot_verbosity):
    dset = IIRrational_data(
        "simple1",
        instance_num=1,
        set_num=1,
    )

    v2_main(["-h"])
    return


def test_prog(request, browser, plotsections, plot_verbosity):
    dset = IIRrational_data(
        "simple1",
        instance_num=1,
        set_num=1,
    )

    print(v2_main(["-z", "+1-1i", "+1+1J"]))
    return
