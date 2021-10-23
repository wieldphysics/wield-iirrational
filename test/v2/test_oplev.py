# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest
from os import path

from IIRrational import v2
from IIRrational.testing.plots import plot_on_assert

from IIRrational.testing.utilities import (
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
def test_oplev_iy(request, browser, plotsections, plot_verbosity):
    from IIRrational.file_io import load_any

    fname = matlab_test_data.matfiles["OplevPlant.mat"]
    fdict = load_any(fname=fname, ftype="mat")

    out = v2.data2filter(
        F_Hz=fdict["ap"]["iy"]["ff"],
        data=fdict["ap"]["iy"]["plant"],
        F_nyquist_Hz=None,
        delay_s=0,
        hints=[sign_validate_and_plot_hint(__file__, request)],
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway=True):
        pass
    return


@skipdeco
def test_oplev_iy(request, browser, plotsections, plot_verbosity):
    from IIRrational.file_io import load_any

    fname = matlab_test_data.matfiles["OplevPlant.mat"]
    fdict = load_any(fname=fname, ftype="mat")

    out = v2.data2filter(
        F_Hz=fdict["ap"]["iy"]["ff"],
        data=fdict["ap"]["iy"]["plant"],
        never_unstable_poles=True,
        never_unstable_zeros=True,
        mode="rational",
        hints=[stability_validate_and_plot_hint(__file__, request)],
    )
    print(out.fitter.poles)
    with plot_on_assert(__file__, request, out.fitter, plot_anyway=True):
        pass
    return
