# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os
import os.path as path
import contextlib

from IIRrational import v1
from IIRrational import fitters_rational, fitters_ZPK
from IIRrational.testing.plots import plot_on_assert

try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False

@pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
def test_quad_M0_P2V_10Hz(request, browser, plotsections, plot_verbosity):
    data_name = 'quad_M0_P2V_10Hz'
    out = v1.data2filter(
        IIRrational_data_dev(
            data_name,
            instance_num = 1,
            set_num = 2,
        ),
        order_initial = 80,
    )

    with plot_on_assert(__file__, request, out.fitter, plot_anyway = True):
        pass
        #assert(False)
    return

@pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
def test_quad_M0_P2V_10Hz_S(request, browser, plotsections, plot_verbosity):
    data_name = 'quad_M0_P2V_10Hz'
    out = v1.data2filter(
        IIRrational_data_dev(
            data_name,
            instance_num = 1,
            set_num = 2,
        ),
        F_nyquist_Hz = None,
        order_initial = 80,
        delay_s = 0,
    )
    with plot_on_assert(__file__, request, out.fitter, plot_anyway = True):
        pass
        #assert(False)
    return

#@pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
#def test_quad_M0_P2V_10Hz_rat_Z(browser, plotsections, plot_verbosity):
#    data_name = 'quad_M0_P2V_10Hz'
#    out = v1.rational_disc_fit(
#        IIRrational_data_dev(
#            data_name,
#            instance_num = 1,
#            set_num = 2,
#        ),
#        nyquist_final_Hz = 1e3,
#        order = 100,
#    )
#
#    import matplotlib.pyplot as plt
#    from IIRrational import plots
#    axB = plots.plot_fitter_flag(
#        out,
#        #don't use jpg as matplotlib can't write it out of the box (needs pillow)
#        fname = path.join(path.split(__file__)[0], 'out-cplx/test-rational-Z.png'),
#    )
#    plt.close(axB.fig)
#    plt.show()
#    return
#
#@pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
#def test_quad_M0_P2V_10Hz_rat_S(browser, plotsections, plot_verbosity):
#    data_name = 'quad_M0_P2V_10Hz'
#    out = v1.rational_disc_fit(
#        IIRrational_data_dev(
#            data_name,
#            instance_num = 1,
#            set_num = 2,
#        ),
#        nyquist_final_Hz = 0,
#        order = 100,
#    )
#
#    import matplotlib.pyplot as plt
#    from IIRrational import plots
#    axB = plots.plot_fitter_flag(
#        out,
#        #don't use jpg as matplotlib can't write it out of the box (needs pillow)
#        fname = path.join(path.split(__file__)[0], 'out-cplx/test-rational-S.png'),
#    )
#    plt.close(axB.fig)
#    plt.show()
#    return


if __name__=='__main__':
    test_simple1(False, None, 10)
