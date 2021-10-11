# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from IIRrational import v2
from IIRrational.testing.plots import plot_on_assert, digest_on_assert

from IIRrational.testing.utilities import (
    sign_validate_and_plot_hint,
)

try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


@pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
def test_HTTS_P_L(request, browser, plotsections, plot_verbosity):
    data = IIRrational_data_dev('HTTS_P_L')
    out = v2.data2filter(
        F_Hz = data.F_Hz,
        data = data.data,
        SNR = data.SNR,
        order_initial = 20,
        hints = [sign_validate_and_plot_hint(__file__, request)],
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway = True):
        pass
    #with digest_on_assert(__file__, request, out, plot_anyway = True):
    #    pass
    return

@pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
def test_HTTS_Y_L(request, browser, plotsections, plot_verbosity):
    data = IIRrational_data_dev('HTTS_Y_L')
    out = v2.data2filter(
        F_Hz = data.F_Hz,
        data = data.data,
        SNR = data.SNR,
        order_initial = 20,
        hints = [sign_validate_and_plot_hint(__file__, request)],
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway = True):
        pass
    return

@pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
def test_HTTS_Y_L_reldeg(request, browser, plotsections, plot_verbosity):
    def reldeg_validate(aid, fitter):
        rdmax = aid.hint('relative_degree_max')
        rdmin = aid.hint('relative_degree_min')
        orders = aid.fitter_orders(fitter)
        #print("DEGREE", orders.total, orders.reldeg, rdmin, rdmax)
        if rdmax is not None:
            assert(orders.reldeg <= rdmax)
        if rdmin is not None:
            assert(orders.reldeg >= rdmin)

    hint = {
        'fitter_update_validate'   : reldeg_validate,
        'fitter_check_validate'    : reldeg_validate,
    }
    data = IIRrational_data_dev('HTTS_Y_L')
    out = v2.data2filter(
        F_Hz = data.F_Hz,
        data = data.data,
        SNR = data.SNR,
        order_initial = 20,
        hints = [hint],
        relative_degree = -2,
        log_level = 10,
        log_level_debug = 10,
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway = True):
        pass
    return

