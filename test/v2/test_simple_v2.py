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
    data_name = "simple1"

    dset = IIRrational_data(
        data_name,
        instance_num=1,
        set_num=1,
    )
    out = v2.data2filter(
        F_Hz=dset.F_Hz,
        SNR=dset.SNR,
        data=dset.data,
        order_initial=10,
        delay_s=0,
        F_nyquist_Hz=None,
        hints=[
            testing.validate_plot_log(__file__, request),
        ],
    )

    folder = path.join(path.split(__file__)[0], "plots")
    axB = out.investigate_order_plot()
    axB.save(path.join(folder, data_name))
    return


def test_simple1_fit(request, browser, plotsections, plot_verbosity):
    data_name = "simple1"
    dset = IIRrational_data(
        data_name,
        instance_num=1,
        set_num=1,
    )
    out = v2.data2filter(
        F_Hz=dset.F_Hz,
        SNR=dset.SNR,
        data=dset.data,
        ZPK=dset.bestfit_ZPK_s,
        order_initial=10,
        delay_s=0,
        F_nyquist_Hz=None,
        mode="fit",
        hints=[
            testing.validate_plot_log(__file__, request),
        ],
    )

    with testing.plot_on_assert(__file__, request, out.fitter, plot_anyway=False):
        pass
    print(out.as_foton_str_ZPKsf())
    return
